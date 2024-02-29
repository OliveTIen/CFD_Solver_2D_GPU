#include "GPUSolver2.h"
#include "../FVM_2D.h"
#include "space/CalculateGradient2.cuh"
#include "space/CalculateFlux2.cuh"
#include "math/PhysicalConvertKernel.h"
#include "time/Evolve.h"
#include "time/CalculateDt.h"
#include "../output/LogWriter.h"

void GPU::GPUSolver2::initialze() {

	// 设置GPU
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

	// 申请内存
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

	// 申请device内存
	try {
		node_device.cuda_alloc(num_node);
		element_device.cuda_alloc(num_element);
		elementField_device.cuda_alloc(num_element);
		edge_device.cuda_alloc(num_edge);
	}
	catch (const char* e) {
		// ! 异常处理未完成
		std::cout << e << std::endl;
		cudaError_t error = cudaError_t::cudaErrorDeviceUninitialized;
		exit(-1);
	}

	// 初始化Element_2D vector中的GPUindex
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
	// 初始化host
	initialize_nodeHost(pFVM2D, num_node);
	initialize_elementHost(pFVM2D, num_element);
	initialize_edgeHost(pFVM2D, num_edge);
	initialize_boundary(pFVM2D);

	//std::cout << "cuda_memcpy\n";
	// 复制host数据到device
	GPU::NodeSoA::cuda_memcpy(&node_device, &node_host, cudaMemcpyHostToDevice);
	GPU::ElementSoA::cuda_memcpy(&element_device, &element_host, cudaMemcpyHostToDevice);
	GPU::FieldSoA::cuda_memcpy(&elementField_device, &elementField_host, cudaMemcpyHostToDevice);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);

}

void GPU::GPUSolver2::iteration(REAL& t, REAL T) {
	// 更新时间t
	update_ruvp_Uold();
	REAL dt = calculateDt(t, T);
	t += dt;
	iterationDevice(dt);
	device_to_host();
	iterationHost(dt);
	host_to_device();
}

void GPU::GPUSolver2::iterationDevice(REAL dt) {


	// 计算单元梯度 对device操作
	GPU::calculateGradient2(this->element_device, this->elementField_device);
	cudaDeviceSynchronize();

	//// 时间推进 对device操作
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
	// 初始化node_host
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
	假设是三角形
	*/

	// 初始化element_host和elementAdjacent
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
	const int nodePerElement = 3;
	const int nValue = 4;
//#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {

		Element_2D& element_i = pFVM2D->elements[i];
		// ID xy volume
		// 注意volume是用的以前的函数计算，如果要重新读取的话，xy volume需要放在后面，因为
		// node尚未全部初始化
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
	// 初始化edge_host
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

		// 计算边界坐标，等于节点坐标平均值
		int tmpNodeID0 = edge_host.nodes[0][i];// node ID
		int tmpNodeID1 = edge_host.nodes[1][i];
		edge_host.xy[0][i] = (node_host.xy[0][tmpNodeID0] + node_host.xy[0][tmpNodeID1]) / 2.0;
		edge_host.xy[1][i] = (node_host.xy[1][tmpNodeID0] + node_host.xy[1][tmpNodeID1]) / 2.0;
		// 计算normal
		//edge_i.getDirectionN(edge_host.normal[0][i], edge_host.normal[1][i]);
		REAL dx = node_host.xy[0][tmpNodeID1] - node_host.xy[0][tmpNodeID0];
		REAL dy = node_host.xy[1][tmpNodeID1] - node_host.xy[1][tmpNodeID0];
		edge_host.normal[0][i] = dy / edge_i.length;
		edge_host.normal[1][i] = -dx / edge_i.length;
	}

}

void GPU::GPUSolver2::initialize_boundary(void* _pFVM2D_) {
	// 遍历边界，找到周期边界，然后遍历它的边，更新elementR
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
	注：经过测试，发现auto指针对ij指标的访问是可行的
	无论是auto还是auto&，均可以
	对于auto，虽然原数组是REAL*a[4]类型，新指针是REAL**类型，但不影响访问
	对于auto&，毋庸置疑可以访问
	*/
	auto& U = this->elementField_host.U;
	auto& U_old = this->element_U_old;
	auto& ruvp = this->element_vruvp;
	for (int i = 0; i < this->elementField_host.num; i++) {
		// 更新host端的U_old
		for (int j = 0; j < 4; j++) {
			U_old[j][i] = U[j][i];
		}
		// 更新host端的ruvp 
		// U:rho,rho_u,rho_v,rho_E ruvp:rho,u,v,p
		ruvp[0][i] = U[0][i];
		ruvp[1][i] = U[1][i] / U[0][i];
		ruvp[2][i] = U[2][i] / U[0][i];
		ruvp[3][i] = ruvp[0][i] * (gamma - 1) * (U[3][i] / U[0][i] - (ruvp[1][i] * ruvp[1][i] + ruvp[2][i] * ruvp[2][i]) * 0.5);
	}
}

REAL GPU::GPUSolver2::calculateDt(REAL t, REAL T) {
	// 目前默认是无粘的，因此ifViscous和ifConstViscous都为false
	return GPU::calculateDt(
		t, T, Constant::gamma, Constant::Re, Constant::Pr, GlobalPara::time::CFL, Constant::R,
		false, false, *this
	);
}

void GPU::GPUSolver2::device_to_host() {
	// 将device端的数据复制到host端
	GPU::NodeSoA::cuda_memcpy(&node_host, &node_device, cudaMemcpyDeviceToHost);
	GPU::ElementSoA::cuda_memcpy(&element_host, &element_device, cudaMemcpyDeviceToHost);
	GPU::FieldSoA::cuda_memcpy(&elementField_host, &elementField_device, cudaMemcpyDeviceToHost);
	GPU::EdgeSoA::cuda_memcpy(&edge_host, &edge_device, cudaMemcpyDeviceToHost);

}

void GPU::GPUSolver2::iterationHost(REAL dt) {
	// 主机端的迭代
	//GPU::calculateFluxHost(element_host, elementField_host, edge_host);
	GPU::Time::EvolveHost(dt, GlobalPara::time::time_advance, element_host, elementField_host, edge_host, edge_periodic_pair);
}

void GPU::GPUSolver2::host_to_device() {
	// 将host端的数据复制到device端
	GPU::NodeSoA::cuda_memcpy(&node_device, &node_host, cudaMemcpyHostToDevice);
	GPU::ElementSoA::cuda_memcpy(&element_device, &element_host, cudaMemcpyHostToDevice);
	GPU::FieldSoA::cuda_memcpy(&elementField_device, &elementField_host, cudaMemcpyHostToDevice);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);
}

void GPU::GPUSolver2::updateOutputNodeField() {

	/*
	可能的隐患：
	这里nNodePerElement为4，
	对于三角形单元，有两种处理方式，第一种让nodeID[3]=-1，另一种让nodeID[3]=nodeID[2]
	当nodeID[3]=-1时，会跳过该节点，防止越界访问数组
	当nodeID[3]=nodeID[2]时，邻居数量会加1，同时值会加两遍。这导致nodeID[2]权重偏大。
	*/

	auto& node_ruvp = this->outputNodeField.ruvp;
	int num_node = this->outputNodeField.num_node;
	int num_element = this->element_host.num_element;
	const int nNodePerElement = 4;
	const int nValuePerNode = 4;
	int* node_neighborElement_num;// 记录每个节点的邻居单元数量
	// 申请资源 初始化
	node_neighborElement_num = new int[num_node]();// 小括号自动初始化为0
	for (int i = 0; i < 4; i++) {
		memset(node_ruvp[i], 0, num_node * sizeof(REAL));
	}

	// 将单元值加到其子节点值上
	for (int iElement = 0; iElement < num_element; iElement++) {
		// 统计每个节点的邻居单元数量，并修改节点值
		for (int jNode = 0; jNode < nNodePerElement; jNode++) {
			// 获取单元节点ID
			int GPUID_of_node = this->element_host.nodes[jNode][iElement];// GPUID of node 0
			// nodeID[3]=-1时，跳过该节点
			if (GPUID_of_node < 0 || GPUID_of_node >= num_node)continue;// 跳过本次循环，但不跳出循环体
			// 该ID对应的邻居单元数量+1
			node_neighborElement_num[GPUID_of_node]++;
			// 该ID对应的节点所有值加上单元值
			for (int kValue = 0; kValue < nValuePerNode; kValue++) {
				node_ruvp[kValue][GPUID_of_node] += this->element_vruvp[kValue][iElement];
			}
		}
	}

	// 节点值除以邻居单元数，得到平均值，作为节点ruvp值
	for (int iNode = 0; iNode < num_node; iNode++) {
		// 为了避免除以0，分母等于0则跳过
		if (node_neighborElement_num[iNode] == 0)continue;
		// node_ruvp除以邻居单元数，得到平均值
		for (int kValue = 0; kValue < nValuePerNode; kValue++) {
			node_ruvp[kValue][iNode] /= node_neighborElement_num[iNode];
		}
	}

	// 释放资源
	delete[] node_neighborElement_num;
}
