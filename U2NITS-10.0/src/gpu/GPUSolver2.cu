#include "GPUSolver2.h"
#include "../FVM_2D.h"
#include "space/CalculateGradient2.cuh"
#include "space/CalculateFlux2.cuh"
#include "math/PhysicalConvertKernel.h"
#include "time/Evolve.h"
#include "time/CalculateDt.h"
#include "../output/LogWriter.h"

void GPU::GPUSolver2::initialze() {

	// ����GPU
	int iDeviceCount = 0;
	cudaError_t error = cudaGetDeviceCount(&iDeviceCount);
	if (error != cudaSuccess || iDeviceCount <= 0) {
		LogWriter::writeLog(cudaGetErrorString(error), LogWriter::Error);
		exit(-1);
	}
	printf("Num of GPU: %d\n", iDeviceCount);

	int iDev = 0;
	error = cudaSetDevice(iDev);
	if (error != cudaSuccess) {
		LogWriter::writeLog(cudaGetErrorString(error), LogWriter::Error);
		exit(-1);
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
	//periodicBoundary.alloc(num_edge);

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
		exit(-1);
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

	//std::cout << "initialize_nodeHost\n";
	// ��ʼ��host
	initialize_nodeHost(pFVM2D, num_node);
	initialize_elementHost(pFVM2D, num_element);
	initialize_edgeHost(pFVM2D, num_edge);
	initialize_boundary(pFVM2D);

	//std::cout << "cuda_memcpy\n";
	// ����host���ݵ�device
	GPU::NodeSoA::cuda_memcpy(&node_device, &node_host, cudaMemcpyHostToDevice);
	GPU::ElementSoA::cuda_memcpy(&element_device, &element_host, cudaMemcpyHostToDevice);
	GPU::FieldSoA::cuda_memcpy(&elementField_device, &elementField_host, cudaMemcpyHostToDevice);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);

}

void GPU::GPUSolver2::iteration(REAL& t, REAL T) {
	// ����ʱ��t
	update_ruvp_Uold();
	REAL dt = calculateDt(t, T);
	t += dt;
	iterationDevice(dt);
	device_to_host();
	iterationHost(dt);
	host_to_device();
}

void GPU::GPUSolver2::iterationDevice(REAL dt) {


	// ���㵥Ԫ�ݶ� ��device����
	GPU::calculateGradient2(this->element_device, this->elementField_device);
	cudaDeviceSynchronize();

	//// ʱ���ƽ� ��device����
	// TODO
	//cudaDeviceSynchronize();
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
	//this->periodicBoundary.free();
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
	/*
	������������
	*/

	// ��ʼ��element_host��elementAdjacent
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
	const int nodePerElement = 3;
	const int nValue = 4;
//#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {

		Element_2D& element_i = pFVM2D->elements[i];
		// ID xy volume
		// ע��volume���õ���ǰ�ĺ������㣬���Ҫ���¶�ȡ�Ļ���xy volume��Ҫ���ں��棬��Ϊ
		// node��δȫ����ʼ��
		element_host.ID[i] = element_i.GPUID;
		element_host.xy[0][i] = element_i.x;
		element_host.xy[1][i] = element_i.y;
		element_host.volume[i] = element_i.calArea(pFVM2D);
		for (int j = 0; j < nodePerElement; j++) {
			element_host.nodes[j][i] = pFVM2D->getNodeByID(element_i.nodes[j])->GPUID;
			element_host.edges[j][i] = element_i.pEdges[j]->GPUID;
		}
		for (int j = 0; j < nValue; j++) {
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
		if (edge_i.pElement_R == nullptr) {
			edge_host.elementR[i] = -1;
		}
		else {
			edge_host.elementR[i] = edge_i.pElement_R->GPUID;
		}
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

void GPU::GPUSolver2::initialize_boundary(void* _pFVM2D_) {
	// �����߽磬�ҵ����ڱ߽磬Ȼ��������ıߣ�����elementR
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
	BoundaryManager_2D& boundaryManager = pFVM2D->boundaryManager;
	for (VirtualBoundarySet_2D& boundary : boundaryManager.boundaries) {
		int bType = boundary.type;
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {
			for (Edge_2D* pEdge : boundary.pEdges) {
				int edgeID = pEdge->GPUID;
				int edgeID_pair = boundaryManager.get_pairEdge_periodic(pEdge)->GPUID;
				int elementL_pair = this->edge_host.elementL[edgeID_pair];
				this->edge_host.elementR[edgeID] = elementL_pair;
				edge_periodic_pair.insert(std::pair<int, int>(edgeID, edgeID_pair));
			}
		}
	}
}

void GPU::GPUSolver2::update_ruvp_Uold() {

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
}

REAL GPU::GPUSolver2::calculateDt(REAL t, REAL T) {
	// ĿǰĬ������ճ�ģ����ifViscous��ifConstViscous��Ϊfalse
	return GPU::calculateDt(
		t, T, Constant::gamma, Constant::Re, Constant::Pr, GlobalPara::time::CFL, Constant::R,
		false, false, *this
	);
}

void GPU::GPUSolver2::device_to_host() {
	// ��device�˵����ݸ��Ƶ�host��
	GPU::NodeSoA::cuda_memcpy(&node_host, &node_device, cudaMemcpyDeviceToHost);
	GPU::ElementSoA::cuda_memcpy(&element_host, &element_device, cudaMemcpyDeviceToHost);
	GPU::FieldSoA::cuda_memcpy(&elementField_host, &elementField_device, cudaMemcpyDeviceToHost);
	GPU::EdgeSoA::cuda_memcpy(&edge_host, &edge_device, cudaMemcpyDeviceToHost);

}

void GPU::GPUSolver2::iterationHost(REAL dt) {
	// �����˵ĵ���
	//GPU::calculateFluxHost(element_host, elementField_host, edge_host);
	GPU::Time::EvolveHost(dt, GlobalPara::time::time_advance, element_host, elementField_host, edge_host, edge_periodic_pair);
}

void GPU::GPUSolver2::host_to_device() {
	// ��host�˵����ݸ��Ƶ�device��
	GPU::NodeSoA::cuda_memcpy(&node_device, &node_host, cudaMemcpyHostToDevice);
	GPU::ElementSoA::cuda_memcpy(&element_device, &element_host, cudaMemcpyHostToDevice);
	GPU::FieldSoA::cuda_memcpy(&elementField_device, &elementField_host, cudaMemcpyHostToDevice);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);
}

void GPU::GPUSolver2::updateOutputNodeField() {

	/*
	���ܵ�������
	����nNodePerElementΪ4��
	���������ε�Ԫ�������ִ���ʽ����һ����nodeID[3]=-1����һ����nodeID[3]=nodeID[2]
	��nodeID[3]=-1ʱ���������ýڵ㣬��ֹԽ���������
	��nodeID[3]=nodeID[2]ʱ���ھ��������1��ͬʱֵ������顣�⵼��nodeID[2]Ȩ��ƫ��
	*/

	auto& node_ruvp = this->outputNodeField.ruvp;
	int num_node = this->outputNodeField.num_node;
	int num_element = this->element_host.num_element;
	const int nNodePerElement = 4;
	const int nValuePerNode = 4;
	int* node_neighborElement_num;// ��¼ÿ���ڵ���ھӵ�Ԫ����
	// ������Դ ��ʼ��
	node_neighborElement_num = new int[num_node]();// С�����Զ���ʼ��Ϊ0
	for (int i = 0; i < 4; i++) {
		memset(node_ruvp[i], 0, num_node * sizeof(REAL));
	}

	// ����Ԫֵ�ӵ����ӽڵ�ֵ��
	for (int iElement = 0; iElement < num_element; iElement++) {
		// ͳ��ÿ���ڵ���ھӵ�Ԫ���������޸Ľڵ�ֵ
		for (int jNode = 0; jNode < nNodePerElement; jNode++) {
			// ��ȡ��Ԫ�ڵ�ID
			int GPUID_of_node = this->element_host.nodes[jNode][iElement];// GPUID of node 0
			// nodeID[3]=-1ʱ�������ýڵ�
			if (GPUID_of_node < 0 || GPUID_of_node >= num_node)continue;// ��������ѭ������������ѭ����
			// ��ID��Ӧ���ھӵ�Ԫ����+1
			node_neighborElement_num[GPUID_of_node]++;
			// ��ID��Ӧ�Ľڵ�����ֵ���ϵ�Ԫֵ
			for (int kValue = 0; kValue < nValuePerNode; kValue++) {
				node_ruvp[kValue][GPUID_of_node] += this->element_vruvp[kValue][iElement];
			}
		}
	}

	// �ڵ�ֵ�����ھӵ�Ԫ�����õ�ƽ��ֵ����Ϊ�ڵ�ruvpֵ
	for (int iNode = 0; iNode < num_node; iNode++) {
		// Ϊ�˱������0����ĸ����0������
		if (node_neighborElement_num[iNode] == 0)continue;
		// node_ruvp�����ھӵ�Ԫ�����õ�ƽ��ֵ
		for (int kValue = 0; kValue < nValuePerNode; kValue++) {
			node_ruvp[kValue][iNode] /= node_neighborElement_num[iNode];
		}
	}

	// �ͷ���Դ
	delete[] node_neighborElement_num;
}
