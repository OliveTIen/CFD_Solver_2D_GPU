#include "GPUSolver2.h"
#include "../FVM_2D.h"
#include "space/CalculateGradient2.h"
#include "space/CalculateFlux2.h"
#include "math/PhysicalConvertKernel.h"

void GPU::GPUSolver2::initialze() {

	// ����GPU
	int iDeviceCount = 0;
	cudaError_t error = cudaGetDeviceCount(&iDeviceCount);
	if (error != cudaSuccess || iDeviceCount <= 0) {
		printf("No GPU found\n");
		throw error;
	}
	printf("Num of GPU: %d\n", iDeviceCount);

	int iDev = 0;
	error = cudaSetDevice(iDev);
	if (error != cudaSuccess) {
		printf("Fail to set GPU %d\n", iDev);
		throw error;
	}
	printf("Activate GPU: %d\n", iDev);

	// �����ڴ�
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;
	const int num_node = (int)pFVM2D->nodes.size();
	const int num_element = (int)pFVM2D->elements.size();
	//const int num_neighbor = 3;
	const int num_edge = (int)pFVM2D->edges.size();
	this->node_host.alloc(num_node);
	this->element_host.alloc(num_element);
	this->elementField_host.alloc(num_element);
	this->edge_host.alloc(num_edge);
	for (int j = 0; j < 4; j++) {
		element_vruvp[j] = new REAL[num_element];
		element_U_old[j] = new REAL[num_element];
	}
	outputNodeField.alloc(num_node);

	// ����device�ڴ�
	try {
		node_device.cuda_alloc(num_node);
		element_device.cuda_alloc(num_element);
		elementField_device.cuda_alloc(num_element);
		edge_device.cuda_alloc(num_edge);
	}
	catch (const char* e) {
		// ! �쳣����δ���
		std::cout << e << std::endl;
		cudaError_t error = cudaError_t::cudaErrorDeviceUninitialized;
		throw error;
	}

	// ��ʼ��Element_2D vector�е�GPUindex
//#pragma omp parallel for
	for (int i = 0; i < num_node; i++) {
		pFVM2D->nodes[i].GPUID = i;
	}
//#pragma omp barrier
//#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {
		pFVM2D->elements[i].GPUID = i;
	}
//#pragma omp barrier
//#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		pFVM2D->edges[i].GPUID = i;
	}
//#pragma omp barrier

	// ��ʼ��host
	initialize_nodeHost(pFVM2D, num_node);
	initialize_elementHost(pFVM2D, num_element);
	initialize_edgeHost(pFVM2D, num_edge);

	// ����host���ݵ�device
	GPU::NodeSoA::cuda_memcpy(&node_device, &node_host, cudaMemcpyHostToDevice);
	GPU::ElementSoA::cuda_memcpy(&element_device, &element_host, cudaMemcpyHostToDevice);
	GPU::FieldSoA::cuda_memcpy(&elementField_device, &elementField_host, cudaMemcpyHostToDevice);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);

}

void GPU::GPUSolver2::iteration() {


	// ���㵥Ԫ�ݶ� ��device����
	GPU::calculateGradient2(this->element_device, this->elementField_device);
	cudaDeviceSynchronize();

	// ����߽���ֵͨ�� ��device����
	GPU::calculateFlux2(this->element_device, this->elementField_device, this->edge_device);
	cudaDeviceSynchronize();
}

void GPU::GPUSolver2::finalize() {
	this->node_host.free();
	this->node_device.cuda_free();
	this->element_host.free();
	this->element_device.cuda_free();
	this->elementField_host.free();
	this->elementField_device.cuda_free();
	this->edge_host.free();
	this->edge_device.cuda_free();
	for (int j = 0; j < 4; j++) {
		delete[] element_vruvp[j];
		delete[] element_U_old[j];
	}
	this->outputNodeField.free();
}

void GPU::GPUSolver2::initialize_nodeHost(void* _pFVM2D_, int num_node) {
	// ��ʼ��node_host
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
//#pragma omp parallel for
	for (int i = 0; i < num_node; i++) {
		Node_2D node_i = pFVM2D->nodes[i];
		node_host.ID[i] = node_i.GPUID;
		node_host.xy[0][i] = node_i.x;
		node_host.xy[1][i] = node_i.y;
	}
}

void GPU::GPUSolver2::initialize_elementHost(void* _pFVM2D_, int num_element) {
	// ��ʼ��element_host��elementAdjacent
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
//#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {

		Element_2D& element_i = pFVM2D->elements[i];

		/*
		���ݽṹ�ο���
			class ElementSoA {
	public:
		int num_element;

		int* ID;
		int* nodes[4];
		int* edges[4];
		int* neighbors[4];
		REAL* x;
		REAL* y;
		REAL* U[4];
		REAL* Ux[4];
		REAL* Uy[4];
		REAL* Flux[4];
		
		*/
		// ID xy volume
		// ע��volume���õ���ǰ�ĺ������㣬���Ҫ���¶�ȡ�Ļ���xy volume��Ҫ���ں��棬��Ϊ
		// node��δȫ����ʼ��
		element_host.ID[i] = element_i.GPUID;
		element_host.xy[0][i] = element_i.x;
		element_host.xy[1][i] = element_i.y;
		element_host.volume[i] = element_i.calArea(pFVM2D);
		for (int j = 0; j < 4; j++) {
			element_host.nodes[j][i] = pFVM2D->getNodeByID(element_i.nodes[j])->GPUID;
			element_host.edges[j][i] = element_i.pEdges[j]->GPUID;
			// ��
			elementField_host.U[j][i] = element_i.U[j];
			elementField_host.Ux[j][i] = element_i.Ux[j];
			elementField_host.Uy[j][i] = element_i.Uy[j];
			elementField_host.Flux[j][i] = element_i.Flux[j];
		}

		std::vector<Element_2D*> neighbors_element_i = element_i.findNeighbor();
		size_t num_neighbor = neighbors_element_i.size();
		for (int j = 0; j < num_neighbor; j++) {
			if (neighbors_element_i[j] == nullptr) {
				element_host.neighbors[j][i] = -1;
			}
			else {
				element_host.neighbors[j][i] = neighbors_element_i[j]->GPUID;
			}
		}
	}


}

void GPU::GPUSolver2::initialize_edgeHost(void* _pFVM2D_, int num_edge) {
	// ��ʼ��edge_host
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
//#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		Edge_2D& edge_i = pFVM2D->edges[i];
		edge_host.ID[i] = edge_i.GPUID;
		edge_host.nodes[0][i] = pFVM2D->getNodeByID(edge_i.nodes[0])->GPUID;
		edge_host.nodes[1][i] = pFVM2D->getNodeByID(edge_i.nodes[1])->GPUID;
		edge_host.setID[i] = edge_i.setID;
		edge_host.elementL[i] = edge_i.pElement_L->GPUID;
		edge_host.elementR[i] = edge_i.pElement_R->GPUID;
		edge_host.length[i] = edge_i.length;
		edge_host.distanceOfElements[i] = edge_i.refLength;

		// ����߽����꣬���ڽڵ�����ƽ��ֵ
		int tmpNodeID0 = edge_host.nodes[0][i];// node ID
		int tmpNodeID1 = edge_host.nodes[1][i];
		edge_host.xy[0][i] = (node_host.xy[0][tmpNodeID0] + node_host.xy[0][tmpNodeID1]) / 2.0;
		edge_host.xy[1][i] = (node_host.xy[1][tmpNodeID0] + node_host.xy[1][tmpNodeID1]) / 2.0;
		// ����normal
		//edge_i.getDirectionN(edge_host.normal[0][i], edge_host.normal[1][i]);
		REAL dx = node_host.xy[0][tmpNodeID1] - node_host.xy[0][tmpNodeID0];
		REAL dy = node_host.xy[1][tmpNodeID1] - node_host.xy[1][tmpNodeID0];
		edge_host.normal[0][i] = dy / edge_i.length;
		edge_host.normal[1][i] = -dx / edge_i.length;
	}
}

void GPU::GPUSolver2::update_host_data(void* _pFVM2D_) {

	REAL gamma = Constant::gamma;
	/*
	ע���������ԣ�����autoָ���ijָ��ķ����ǿ��е�
	������auto����auto&��������
	����auto����Ȼԭ������REAL*a[4]���ͣ���ָ����REAL**���ͣ�����Ӱ�����
	����auto&����ӹ���ɿ��Է���
	*/
	auto& U = this->elementField_host.U;
	auto& U_old = this->element_U_old;
	auto& ruvp = this->element_vruvp;
	for (int i = 0; i < this->elementField_host.num; i++) {
		// ����host�˵�U_old
		for (int j = 0; j < 4; j++) {
			U_old[j][i] = U[j][i];
		}
		// ����host�˵�ruvp 
		// U:rho,rho_u,rho_v,rho_E ruvp:rho,u,v,p
		ruvp[0][i] = U[0][i];
		ruvp[1][i] = U[1][i] / U[0][i];
		ruvp[2][i] = U[2][i] / U[0][i];
		ruvp[3][i] = ruvp[0][i] * (gamma - 1) * (U[3][i] / U[0][i] - (ruvp[1][i] * ruvp[1][i] + ruvp[2][i] * ruvp[2][i]) * 0.5);
	}


	/*
	�ο����룺
		ruvp[0] = U[0];
		ruvp[1] = U[1] / U[0];
		ruvp[2] = U[2] / U[0];
		ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
	*/
	/*
	����autoָ�룺
	const int nVar = 4;
	const int nElement = 8;
	double* a[nVar]{};
	for (int i = 0; i < nVar; i++) {
		a[i] = new double[nElement];
	}
	for (int i = 0; i < nVar; i++) {
		for (int j = 0; j < nElement; j++) {
			a[i][j] = i * 4 + j;
		}
	}
	for (int i = 0; i < nVar; i++) {
		for (int j = 0; j < nElement; j++) {
			std::cout << a[i][j] << " \t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	// ����ָ������ Ҳ������std::ref
	auto a_ref = a;
	for (int i = 0; i < nVar; i++) {
		for (int j = 0; j < nElement; j++) {
			a_ref[i][j] *= 0.1;
		}
	}

	for (int i = 0; i < nVar; i++) {
		for (int j = 0; j < nElement; j++) {
			std::cout << a[i][j] << " \t";
		}
		std::cout << std::endl;
	}

	for (int i = 0; i < nVar; i++) {
		delete[] a[i];
	}
	
	*/
}

void GPU::GPUSolver2::device_to_host() {
	// ��device�˵����ݸ��Ƶ�host��
	GPU::NodeSoA::cuda_memcpy(&node_host, &node_device, cudaMemcpyDeviceToHost);
	GPU::ElementSoA::cuda_memcpy(&element_host, &element_device, cudaMemcpyDeviceToHost);
	GPU::FieldSoA::cuda_memcpy(&elementField_host, &elementField_device, cudaMemcpyDeviceToHost);
	GPU::EdgeSoA::cuda_memcpy(&edge_host, &edge_device, cudaMemcpyDeviceToHost);

}
