#include "GPUSolver.h"
#include "../FVM_2D.h"
#include "space/CalculateGradient.h"
#include "space/CalculateFlux.h"
#include "datatype/EdgeSoA.h"

void GPUSolver::initialze() {
	/*
	将FVM_2D数据转化为Element结构并初始化
	初始化的数据有：
	    单元坐标、U、单元邻居坐标、单元邻居U
	暂时未初始化：
		单元节点坐标、边

	读取文件后，将FVM_2D的数据结构转换为Element结构：
	1. 初始化selfxy, node1, node2, node3 （由于遗留问题，self是在CPU上计算的，以后会转移到GPU上）
	2. 初始化selfValue.U1~U4, neighbor1.U1~U4, neighbor2.U1~U4, neighbor3.U1~U4
	3. 梯度初始化为0


	cudaMalloc应放在循环外，以减少开销 参加书P364
	但是计算单元和计算边是两套不同大小的数组，因此需要申请两套
	主机复制到设备

	需要初始化：
		单元坐标
		单元节点坐标

	*/
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
	const int num_element = (int) pFVM2D->elements.size();
	const int num_neighbor = 3;
	const int num_edge = pFVM2D->edges.size();
	this->element_host.alloc(num_element);
	this->elementAdjacent.alloc(num_neighbor, num_element);
	this->edge_host.alloc(num_edge);

	// 申请device内存
	try{ 
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
		pFVM2D->elements[i].GPUindex = i;
	}
	#pragma omp barrier
	#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		pFVM2D->edges[i].GPUID = i;
	}
	#pragma omp barrier
	
	// 初始化element_host和elementAdjacent
	initialize_elementHost_elementAdjacent(pFVM2D, num_element, num_neighbor);
	// 初始化edge_host
	initialize_edgeHost(pFVM2D, num_edge);

	// 复制host数据到device
	element_host.cuda_copy_to_device(&element_device);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);
}

void GPUSolver::iteration() {
	/*
    
	1. 初始化
		单元数值通量赋为0
		根据单元U初始化单元邻居U

    2. 计算单元梯度 - GPU
	    2.1
        根据坐标、U，计算梯度
        限制器修正梯度
        检查异常值
        输入：单元坐标、单元U、单元邻居坐标、单元邻居U
        输出：单元梯度

	    2.2. 单元梯度复制到CPU，然后更新邻居梯度
		
    4. 计算边界数值通量
        计算边法向量nx ny、边坐标edgex, edgey、边长edgeLength
		    需要参数：节点坐标
        计算边左右单元在边界的U值(左右极限UL、UR)
        根据边类型，修正UL、UR
        计算黎曼问题
            需要的参数：UL UR nx ny edgeLength
            输出：边界的数值通量
    5. 计算单元数值通量
        将边界数值通量加到单元数值通量上 - 规约操作，OpenMP
		(注意每次循环需要将单元数值通量初始化为0)
    6. 显式时间推进

 
	*/


	// --- 计算单元梯度 --- 
	// 输入：单元坐标、单元U、单元邻居坐标、单元邻居U
	// 输出：单元梯度
	GPU::calculateGradient(this->element_device);

	// 更新邻居梯度 此处直接对device操作 根据书P312 可以尝试用整数映射
	GPU::updateNeighborGradient(this->element_device, this->elementAdjacent);

	/*
	// 复制device数据到host
	element_device.cuda_copy_to_host(&element_host);
	// 限制器 目前不需要加
	
	// 更新邻居梯度 1.已添加elementAdjacent数组，存储邻居index信息
	GPU::updateNeighborGradient(this->element_host, this->elementAdjacent);

	// 复制梯度数据到device
	for (int i = 0; i < 3; i++) {
		cudaMemcpy(element_device.neighbors[i].Ux1, element_host.neighbors[i].Ux1, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Ux2, element_host.neighbors[i].Ux2, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Ux3, element_host.neighbors[i].Ux3, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Ux4, element_host.neighbors[i].Ux4, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Uy1, element_host.neighbors[i].Uy1, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Uy2, element_host.neighbors[i].Uy2, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Uy3, element_host.neighbors[i].Uy3, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Uy4, element_host.neighbors[i].Uy4, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
	}
	*/


	// --- 计算边界数值通量 --- 
	// 输入：单元
	// 输出：边界数值通量
	GPU::calculateFlux(this->element_device);
}

// inline函数只能在头文件中定义，否则会报错LNK2001
void GPUSolver::finalize() {
	this->element_host.free(); 
	this->element_device.cuda_free();
	this->elementAdjacent.free();
}

void GPUSolver::initialize_elementHost_elementAdjacent(void* _pFVM2D_, int num_element, int num_neighbor) {
	// 初始化element_host和elementAdjacent
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {

		Element_2D& element_i = pFVM2D->elements[i];

		// 初始化单元坐标、单元U、单元梯度
		element_host.self.isNull[i] = false;
		element_host.self.x[i] = element_i.x;
		element_host.self.y[i] = element_i.y;
		element_host.self.U1[i] = element_i.U[0];
		element_host.self.U2[i] = element_i.U[1];
		element_host.self.U3[i] = element_i.U[2];
		element_host.self.U4[i] = element_i.U[3];
		element_host.self.Ux1[i] = 0;
		element_host.self.Uy1[i] = 0;
		element_host.self.Ux2[i] = 0;
		element_host.self.Uy2[i] = 0;
		element_host.self.Ux3[i] = 0;
		element_host.self.Uy3[i] = 0;
		element_host.self.Ux4[i] = 0;
		element_host.self.Uy4[i] = 0;

		// 初始化邻居坐标、邻居U、邻居梯度
		std::vector<Element_2D*> neighbors_element_i = element_i.findNeighbor();// 邻居指针
		// 判断邻居是否为nullptr，并赋值。此处默认为3个邻居 ! 可能有隐患
		for (int j = 0; j < num_neighbor; j++) {
			if (neighbors_element_i[j] == nullptr) {
				element_host.neighbors[j].isNull[i] = true;
				elementAdjacent.neighbors[j].index[i] = -1;
			}
			else {
				element_host.neighbors[j].isNull[i] = false;
				element_host.neighbors[j].x[i] = neighbors_element_i[j]->x;
				element_host.neighbors[j].y[i] = neighbors_element_i[j]->y;
				element_host.neighbors[j].U1[i] = neighbors_element_i[j]->U[0];
				element_host.neighbors[j].U2[i] = neighbors_element_i[j]->U[1];
				element_host.neighbors[j].U3[i] = neighbors_element_i[j]->U[2];
				element_host.neighbors[j].U4[i] = neighbors_element_i[j]->U[3];

				element_host.neighbors[j].Ux1[i] = 0;
				element_host.neighbors[j].Uy1[i] = 0;
				element_host.neighbors[j].Ux2[i] = 0;
				element_host.neighbors[j].Uy2[i] = 0;
				element_host.neighbors[j].Ux3[i] = 0;
				element_host.neighbors[j].Uy3[i] = 0;
				element_host.neighbors[j].Ux4[i] = 0;
				element_host.neighbors[j].Uy4[i] = 0;

				// 初始化elementAdjacent
				elementAdjacent.neighbors[j].index[i] = neighbors_element_i[j]->GPUindex;
			}


		}
	}

}

void GPUSolver::initialize_edgeHost(void* _pFVM2D_, int num_edge) {
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		Edge_2D& edge_i = pFVM2D->edges[i];
		edge_host.ID[i]=edge_i.GPUID;
		edge_host.nodes[0][i] = edge_i.nodes[0];
		edge_host.nodes[1][i] = edge_i.nodes[1];
		edge_host.setID[i] = edge_i.setID;
		edge_host.elementL[i] = edge_i.pElement_L->GPUindex;
		edge_host.elementR[i] = edge_i.pElement_R->GPUindex;
		edge_host.length[i] = edge_i.length;
		edge_host.distanceOfElements[i] = edge_i.refLength;
	}

	/*
		int* ID;
		int* nodes[2];
		int* setID;
		int* elementL;
		int* elementR;
		REAL* length;
		REAL* distanceOfElements;// 两侧单元中心距离
	
	*/
}
