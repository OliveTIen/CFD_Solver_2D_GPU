
#include "GPUSolver2.h"
#include "OldDataConverter.h"
#include <sstream>
#include "../gpu/GPUGlobalFunction.h"
#include "../FVM_2D.h"
#include "../space/gradient/GradientGPU.h"
#include "../space/gradient/Gradient.h"
#include "../space/FluxGPU.h"
#include "../space/restrict/Restrict.h"
#include "../math/Math.h"
#include "../time/CalculateDt.h"
#include "../time/Evolve.h"
#include "../time/EvolveGPU.h"
#include "../output/LogWriter.h"
#include "../output/ResidualCalculator.h"


void GPU::GPUSolver2::initialize() {

	initializeHost();// 申请内存，用FVM_2D旧数据初始化
	if (GlobalPara::basic::useGPU) {
		initializeDevice();
	}

}

void GPU::GPUSolver2::initializeHost() {
	// 申请内存，用FVM_2D旧数据初始化
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	U2NITS::OldDataConverter oldDataConverter(pFVM2D, this);
	const int num_node = (int)pFVM2D->nodes.size();
	const int num_element = (int)pFVM2D->elements.size();
	const int num_edge = (int)pFVM2D->edges.size();
	const int num_boundary = (int)pFVM2D->boundaryManager.boundaries.size();

	allocateHostMemory(num_node, num_element, num_edge, num_boundary);
	oldDataConverter.Convert_FVM2D_to_HostData();
}

void GPU::GPUSolver2::initializeDevice() {
	// 申请设备内存，将主机复制到设备
	if (!hostMemoryAllocated) {
		LogWriter::logAndPrintError("Host not allocated. Cannot initialize device.\n");
		exit(-1);
	}
	if (!GlobalPara::basic::useGPU) {
		LogWriter::logAndPrintError("useGPU = false. Cannot initialize device.\n");
		exit(-1);
	}

	FVM_2D* pFVM2D = FVM_2D::getInstance();
	const int num_node = (int)pFVM2D->nodes.size();
	const int num_element = (int)pFVM2D->elements.size();
	const int num_edge = (int)pFVM2D->edges.size();
	const int num_boundary = (int)pFVM2D->boundaryManager.boundaries.size();

	setGPUDevice();// 设置GPU
	allocateDeviceMemory(num_node, num_element, num_edge, num_boundary);
	host_to_device();// 复制host数据到device

}

void GPU::GPUSolver2::iteration(real& t, real T) {
	int useGPU = GlobalPara::basic::useGPU;
	if (useGPU == 0) {
		iterationHost(t, T);
	}
	else if (useGPU == 1) {
		iterationGPU(t, T);
	}
	else if (useGPU == 2) {
		iterationHalfGPU(t, T);
	}
	else if (useGPU == 3) {
		iterationHost_2_addRK3(t, T);
	}
	else {
		exit(-1);
	}
}

void GPU::GPUSolver2::iterationHost(REAL& t, REAL T) {
	update_ruvp_Uold();

	real dt = U2NITS::Time::calculateGlobalDt(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::time::CFL, GlobalPara::constant::R,
		element_host, edge_host, element_vruvp);

	t += dt;
	U2NITS::Space::Gradient::Gradient(element_host, node_host, edge_host, elementField_host);
	U2NITS::Time::EvolveHost(dt, GlobalPara::time::time_advance, element_host, elementField_host, edge_host, edge_periodic_pair);
	U2NITS::Space::Restrict::modifyElementFieldU2d(element_host, elementField_host);
}

void GPU::GPUSolver2::iterationHost_2_addRK3(REAL& t, REAL T) {
	update_ruvp_Uold();
	U2NITS::Time::EvolveHost_2_addRK3(t, T, element_host, elementField_host, edge_host, node_host, element_vruvp);
	U2NITS::Space::Restrict::modifyElementFieldU2d(element_host, elementField_host);
}

void GPU::GPUSolver2::iterationHalfGPU(REAL& t, REAL T) {
	// 更新时间t
	update_ruvp_Uold();
	real dt = U2NITS::Time::calculateGlobalDt(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::time::CFL, GlobalPara::constant::R,
		element_host, edge_host, element_vruvp);
	t += dt;

	GPU::Space::Gradient::Gradient_2(element_device, elementField_device, edge_device);
	cudaDeviceSynchronize();
	catchCudaErrorAndExit();

	device_to_host();
	U2NITS::Time::EvolveHost(dt, GlobalPara::time::time_advance, element_host, elementField_host, edge_host, edge_periodic_pair);
	U2NITS::Space::Restrict::modifyElementFieldU2d(element_host, elementField_host);
	host_to_device();

}

void GPU::GPUSolver2::iterationGPU(REAL& t, REAL T) {
	// 更新时间t @host
	update_ruvp_Uold();
	real dt = U2NITS::Time::calculateGlobalDt(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::time::CFL, GlobalPara::constant::R,
		element_host, edge_host, element_vruvp);
	t += dt;

	// device
	GPU::Space::Gradient::Gradient_2(element_device, elementField_device, edge_device);
	catchCudaErrorAndExit();// cudaErrorIllegalAddress(700). an illegal memory access was encountered. 可以尝试用更小的块？
	cudaDeviceSynchronize();
	catchCudaErrorAndExit();

	LogWriter::logAndPrintError("Unimplemented.GPU::Time::EvolveDevice @ GPUSolver2::iterationGPU\n");
	exit(-11);
	GPU::Time::EvolveDevice(dt, GlobalPara::time::time_advance, element_device, elementField_device, edge_device, boundary_device, sDevicePara);
	catchCudaErrorAndExit();
	cudaDeviceSynchronize();
	catchCudaErrorAndExit();

	device_to_host();// 传递到主机，用于输出
	//host_to_device();// 不需要，因为已经将主机与设备同步

}

void GPU::GPUSolver2::iterationDevice(REAL dt) {
	// 未使用

	// 计算单元梯度 对device操作
	GPU::Space::Gradient::Gradient_2(element_device, elementField_device, edge_device);
	//GPU::calculateGradient_old(element_device, elementField_device);
	cudaDeviceSynchronize();
	GPU::Time::EvolveDevice(dt, GlobalPara::time::time_advance, element_device, elementField_device, edge_device, boundary_device, sDevicePara);
	//// 时间推进 对device操作
	// TODO
	//cudaDeviceSynchronize();
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
	U2NITS::Time::EvolveHost(dt, GlobalPara::time::time_advance, element_host, elementField_host, edge_host, edge_periodic_pair);
}

void GPU::GPUSolver2::host_to_device() {
	// 将host端的数据复制到device端
	GPU::NodeSoA::cuda_memcpy(&node_device, &node_host, cudaMemcpyHostToDevice);
	GPU::ElementSoA::cuda_memcpy(&element_device, &element_host, cudaMemcpyHostToDevice);
	GPU::FieldSoA::cuda_memcpy(&elementField_device, &elementField_host, cudaMemcpyHostToDevice);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);
}

void GPU::GPUSolver2::freeMemory() {
	freeHostMemory();
	if (GlobalPara::basic::useGPU) {
		freeDeviceMemory();
	}
}

void GPU::GPUSolver2::freeHostMemory() {
	if (!hostMemoryAllocated) {
		LogWriter::logAndPrintError("cannot free host memory.\n");
		exit(-1);
	}
	// host host/device双向数据
	this->node_host.free();
	this->element_host.free();
	this->elementField_host.free();
	this->edge_host.free();
	this->boundary_host.free();
	// host host数据
	for (int j = 0; j < 4; j++) {
		delete[] element_vruvp[j];
		delete[] element_U_old[j];
	}
	this->outputNodeField.free();
	// 更新状态
	hostMemoryAllocated = false;
}

void GPU::GPUSolver2::freeDeviceMemory() {
	if (!deviceMemoryAllocated) {
		LogWriter::logAndPrintError("cannot free device memory.\n");
		exit(-1);
	}
	// device 双向数据
	this->node_device.cuda_free();
	this->element_device.cuda_free();
	this->elementField_device.cuda_free();
	this->edge_device.cuda_free();
	this->boundary_device.cuda_free();

	// device device数据
	delete this->infPara_device;
	// 更新状态
	deviceMemoryAllocated = false;
}

void GPU::GPUSolver2::updateResidual() {
	ResidualCalculator::cal_residual_GPU(element_U_old, elementField_host, ResidualCalculator::NORM_INF, residualVector);
}

void GPU::GPUSolver2::getResidual() {
	ResidualCalculator::get_residual_functionF(elementField_host, residualVector, ResidualCalculator::NORM_INF);
}

void GPU::GPUSolver2::setGPUDevice() {
	int nDevice = 0;
	int iDevice = 0;	
	cudaError_t error;
	// 获取设备数
	error = cudaGetDeviceCount(&nDevice);
	if (error != cudaSuccess || nDevice <= 0) {
		LogWriter::log(cudaGetErrorString(error), LogWriter::Error);
		exit(-1);
	}
	// 设置工作设备
	error = cudaSetDevice(iDevice);
	if (error != cudaSuccess) {
		LogWriter::log(cudaGetErrorString(error), LogWriter::Error);
		exit(-1);
	}
	std::stringstream ss;
	ss << "Num of GPU: " << nDevice << ", " << "active GPU: " << iDevice << "\n";
	LogWriter::logAndPrint(ss.str());
	// 更新状态
	deviceReady = true;
}

void GPU::GPUSolver2::allocateMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
	allocateHostMemory(num_node, num_element, num_edge, num_boundary);
	allocateDeviceMemory(num_node, num_element, num_edge, num_boundary);
}

void GPU::GPUSolver2::allocateHostMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
	// 申请host内存
	this->node_host.alloc(num_node);
	this->element_host.alloc(num_element);
	this->elementField_host.alloc(num_element);
	this->edge_host.alloc(num_edge);
	for (int j = 0; j < 4; j++) {
		element_vruvp[j] = new REAL[num_element];
		element_U_old[j] = new REAL[num_element];
	}
	boundary_host.alloc(num_boundary);
	outputNodeField.alloc(num_node);
	// 更新状态
	hostMemoryAllocated = true;
}

void GPU::GPUSolver2::allocateDeviceMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
	if (!deviceReady) {
		LogWriter::logAndPrintError("device isn't ready. Cannot allocate device memory\n");
		exit(-1);
	}
	// 申请device内存
	try {
		// 双向拷贝数据
		node_device.cuda_alloc(num_node);
		element_device.cuda_alloc(num_element);
		elementField_device.cuda_alloc(num_element);
		edge_device.cuda_alloc(num_edge);
		boundary_device.cuda_alloc(num_boundary);

		// 仅传入GPU的数据
		using namespace GlobalPara::boundaryCondition::_2D;// using有作用域，不用担心混淆
		infPara_device = new GPU::DGlobalPara(inf::ruvp, inlet::ruvp, outlet::ruvp, 4
			, &GlobalPara::constant::R, &GlobalPara::constant::gamma, &GlobalPara::inviscid_flux_method::flux_conservation_scheme);
		sDevicePara.initialize(
			GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::c0, GlobalPara::constant::gamma, GlobalPara::constant::epsilon, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::constant::mu,
			GlobalPara::inviscid_flux_method::flag_reconstruct, GlobalPara::inviscid_flux_method::flag_gradient,
			GlobalPara::boundaryCondition::_2D::inf::ruvp, GlobalPara::boundaryCondition::_2D::inlet::ruvp, GlobalPara::boundaryCondition::_2D::outlet::ruvp,
			GlobalPara::inviscid_flux_method::flux_conservation_scheme, GlobalPara::inviscid_flux_method::flux_limiter
		);
		
		// 更新状态
		deviceMemoryAllocated = true;
	}
	catch (const char* e) {
		// ! 异常处理未完成
		std::cout << e << std::endl;
		cudaError_t error = cudaError_t::cudaErrorDeviceUninitialized;
		exit(-1);
	}
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
		edge_host.periodicPair[i] = -1;
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


	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
	BoundaryManager_2D& boundaryManager = pFVM2D->boundaryManager;
	// 周期边界
	/**
	* 周期边界有多种处理方法
	* 一种是在EdgeSoA新增成员int periodicPair，优点是查找快，缺点是占用空间大。
	* 另一种是用map或hashTable存储，但是数据结构编写比较复杂，且GPU上运行比较麻烦，
	* 因为要重新编写__device__函数用于查找
	*
	* 目前GPU采用第一种方法
	*/
	for (VirtualBoundarySet_2D& boundary : boundaryManager.boundaries) {
		int bType = boundary.type;
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {
			// 若某边界集set是周期边界，遍历其edge
			for (Edge_2D* pEdge : boundary.pEdges) {
				int edgeID = pEdge->GPUID;
				int edgeID_pair = boundaryManager.get_pairEdge_periodic(pEdge)->GPUID;
				int elementL_pair = this->edge_host.elementL[edgeID_pair];
				// 用pair的elementL更新自己的elementR
				this->edge_host.elementR[edgeID] = elementL_pair;
				// 用map存储pair
				edge_periodic_pair.insert(std::pair<int, int>(edgeID, edgeID_pair));
				// 用EdgeSoA成员存储pair
				this->edge_host.periodicPair[edgeID] = edgeID_pair;
			}
		}
	}

	// boundarySet
	for (int i = 0; i < boundaryManager.boundaries.size(); i++) {
		boundary_host.type[i] = boundaryManager.boundaries[i].type;
	}
}

void GPU::GPUSolver2::update_ruvp_Uold() {
	/*
	用U更新element_vruvp.element_vruvp用于计算时间步长Dt以及输出流场
	*/
	REAL gamma = GlobalPara::constant::gamma;
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
			if (isnan(U[j][i])) {
				LogWriter::logAndPrintError("isnan.\n");
			}
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
