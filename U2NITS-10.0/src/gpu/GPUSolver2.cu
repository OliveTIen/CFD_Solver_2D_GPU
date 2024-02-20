#include "GPUSolver2.h"
#include "../FVM_2D.h"
#include "space/CalculateGradient2.h"
#include "space/CalculateFlux2.h"

void GPU::GPUSolver2::initialze() {

	// 设置GPU
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

	// 申请内存
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;
	const int num_element = (int)pFVM2D->elements.size();
	const int num_neighbor = 3;
	const int num_edge = (int)pFVM2D->edges.size();
	this->element_host.alloc(num_element);
	this->edge_host.alloc(num_edge);

	// 申请device内存
	try {
		element_device.cuda_alloc(num_element);
		edge_device.cuda_alloc(num_edge);
	}
	catch (const char* e) {
		// ! 异常处理未完成
		std::cout << e << std::endl;
		cudaError_t error = cudaError_t::cudaErrorDeviceUninitialized;
		throw error;
	}

	// 初始化Element_2D vector中的GPUindex
#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {
		pFVM2D->elements[i].GPUID = i;
	}
#pragma omp barrier
#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		pFVM2D->edges[i].GPUID = i;
	}
#pragma omp barrier

	// 初始化host
	initialize_elementHost(pFVM2D, num_element);
	initialize_edgeHost(pFVM2D, num_edge);

	// 复制host数据到device
	GPU::ElementSoA::cuda_memcpy(&element_device, &element_host, cudaMemcpyHostToDevice);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);

}

void GPU::GPUSolver2::iteration() {


	// --- 计算单元梯度 --- 
	// 输入：单元坐标、单元U、单元邻居坐标、单元邻居U
	// 输出：单元梯度
	GPU::calculateGradient2(this->element_device);

	// --- 计算边界数值通量 --- 
	// 输入：单元
	// 输出：边界数值通量
	
	GPU::calculateFlux2(this->element_device, this->edge_device);
}

void GPU::GPUSolver2::finalize() {
	this->element_host.free();
	this->element_device.cuda_free();
	this->edge_host.free();
	this->edge_device.cuda_free();

}

void GPU::GPUSolver2::initialize_elementHost(void* _pFVM2D_, int num_element) {
	// 初始化element_host和elementAdjacent
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {

		Element_2D& element_i = pFVM2D->elements[i];

		/*
		数据结构参考：
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

		element_host.ID[i] = element_i.GPUID;
		element_host.x[i] = element_i.x;
		element_host.y[i] = element_i.y;
		for (int j = 0; j < 4; j++) {
			element_host.nodes[j][i] = element_i.nodes[j];
			element_host.edges[j][i] = element_i.pEdges[j]->GPUID;
			// neighbors留在后面
			element_host.U[j][i] = element_i.U[j];
			element_host.Ux[j][i] = element_i.Ux[j];
			element_host.Uy[j][i] = element_i.Uy[j];
			element_host.Flux[j][i] = element_i.Flux[j];
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
	// 初始化edge_host
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		Edge_2D& edge_i = pFVM2D->edges[i];
		edge_host.ID[i] = edge_i.GPUID;
		edge_host.nodes[0][i] = edge_i.nodes[0];
		edge_host.nodes[1][i] = edge_i.nodes[1];
		edge_host.setID[i] = edge_i.setID;
		edge_host.elementL[i] = edge_i.pElement_L->GPUID;
		edge_host.elementR[i] = edge_i.pElement_R->GPUID;
		edge_host.length[i] = edge_i.length;
		edge_host.distanceOfElements[i] = edge_i.refLength;
	}
}
