
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
#include "../math/BasicAlgorithmGPU.h"
#include "../math/CommonValue.h"
#include "../time/CalculateDt.h"
#include "../time/Evolve.h"
#include "../time/EvolveGPU.h"
#include "../output/LogWriter.h"
//#include "../output/ResidualCalculator.h"
#include "../global/CExit.h"
#include "../time/CIteration.h"
#include "../output/FieldWriter.h"
#include "../drivers/CDriver.h"
#include "../output/residual/Residual.h"

GPU::GPUSolver2* GPU::GPUSolver2::pInstance = nullptr;

GPU::GPUSolver2::GPUSolver2() {
	if (pInstance == nullptr) {
		pInstance = this;
	}
	else {
		LogWriter::logAndPrintError("only allow one instance to exist. @GPUSolver2()\n");
		exit(-1);
	}
}

GPU::GPUSolver2::~GPUSolver2() {
	if (pInstance != nullptr) {
		pInstance = nullptr;
	}
	else {
		LogWriter::logAndPrintError("pInstance == nullptr on destruction. @~GPUSolver2()\n");
	}
}

GPU::GPUSolver2* GPU::GPUSolver2::getInstance() {
	if (pInstance == nullptr) {
		LogWriter::logAndPrintError("pInstance == nullptr. @GPUSolver2::getInstance()\n");
		exit(-1);
	}
	return pInstance;
}

void GPU::GPUSolver2::allocateMemory() {
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	const int num_node = (int)pFVM2D->nodes.size();
	const int num_element = (int)pFVM2D->elements.size();
	const int num_edge = (int)pFVM2D->edges.size();
	const int num_boundary = (int)pFVM2D->boundaryManager.boundaries.size();
	if (num_node == 0 || num_element == 0 || num_edge == 0 || num_boundary == 0) {
		LogWriter::logAndPrintError("has not read field, cannot determine how much memory to allocate. @GPU::GPUSolver2::allocateMemory()\n");
		exit(-1);
	}
	//allocateMemory_kernel(num_node, num_element, num_edge, num_boundary);

	allocateHostMemory(num_node, num_element, num_edge, num_boundary);
	if (GlobalPara::basic::useGPU) {
		setGPUDevice();
		allocateDeviceMemory(num_node, num_element, num_edge, num_boundary);
	}
}

//void GPU::GPUSolver2::allocateMemory_kernel(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
//	allocateHostMemory(num_node, num_element, num_edge, num_boundary);
//	if (GlobalPara::basic::useGPU) {
//		setGPUDevice();
//		allocateDeviceMemory(num_node, num_element, num_edge, num_boundary);
//	}
//}

void GPU::GPUSolver2::allocateHostMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
	if (hostMemoryAllocated) {
		LogWriter::logAndPrintError("Has allocated host memory.\n");
		CExit::saveAndExit(-1);
	}
	// 申请host内存
	this->node_host.alloc(num_node);
	this->element_host.alloc(num_element);
	this->elementField_host.alloc(num_element);
	this->edge_host.alloc(num_edge);
	for (int j = 0; j < 4; j++) {
		this->element_vruvp[j] = new myfloat[num_element];
		this->element_U_old[j] = new myfloat[num_element];
	}
	this->boundary_host.alloc(num_boundary);
	FieldWriter::getInstance()->allocNodeFieldDataUsingOutputScheme(this->outputNodeField, num_node);
	this->edgeField_host.alloc(num_edge);
	this->elementFieldVariable_dt_host.alloc(num_element);

	// 更新状态
	hostMemoryAllocated = true;
}

void GPU::GPUSolver2::allocateDeviceMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
	if (deviceMemoryAllocated) {
		LogWriter::logAndPrintError("Has allocated device memory，cannot allocate twice. @GPUSolver2::allocateDeviceMemory\n");
		CExit::saveAndExit(-1);
	}


	if (!hasSetGPUDevice) {
		LogWriter::logAndPrintError("Hasn't set device. Cannot allocate device memory. @GPUSolver2::allocateDeviceMemory\n");
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
		//infPara_device = new GPU::DGlobalPara(inf::ruvp, inlet::ruvp, outlet::ruvp, 4
		//	, &GlobalPara::constant::R, &GlobalPara::constant::gamma, &GlobalPara::inviscid_flux_method::flux_conservation_scheme);

		// 注意：目前mu0不是最新的参数
		sDevicePara.initialize(
			GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::c0, GlobalPara::constant::gamma, U2NITS::Math::EPSILON, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::constant::mu0,
			GlobalPara::inviscid_flux_method::flag_reconstruct, GlobalPara::inviscid_flux_method::flag_gradient,
			GlobalPara::boundaryCondition::_2D::inf::ruvp, GlobalPara::boundaryCondition::_2D::inlet::ruvp, GlobalPara::boundaryCondition::_2D::outlet::ruvp,
			GlobalPara::inviscid_flux_method::flux_conservation_scheme, GlobalPara::inviscid_flux_method::flux_limiter
		);

		edgeField_device.cuda_alloc(num_edge);
		elementFieldVariable_dt_device.cuda_alloc(num_element);
		// 将alphaC中为了规约而补齐的元素进行赋值。规约是取最小值，因此补齐的元素取较大的值
		GPU::Math::assign_elements_in_array_device(num_element, elementFieldVariable_dt_device.num_reduce, elementFieldVariable_dt_device.alphaC, U2NITS::Math::BIG_FLOAT_NORMAL);

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

void GPU::GPUSolver2::initializeData_byOldData() {

	initializeHostData_byOldData();// 申请内存，用FVM_2D旧数据初始化
	if (GlobalPara::basic::useGPU) {
		initializeDeviceData_byHostData();
	}

}

void GPU::GPUSolver2::initializeHostData_byOldData() {
	if (!hostMemoryAllocated) {
		LogWriter::logAndPrintError("Host data not allocated. @GPU::GPUSolver2::initializeHostData_byOldData()\n");
		exit(-1);
	}
	if (hostDataInitialized) {
		LogWriter::logAndPrintError("Host data has initialized. @GPU::GPUSolver2::initializeHostData_byOldData()\n");
		exit(-1);
	}

	// 用FVM_2D旧数据初始化host数据
	U2NITS::OldDataConverter oldDataConverter(FVM_2D::getInstance(), this);
	oldDataConverter.Convert_FVM2D_to_HostData();

	hostDataInitialized = true;
}

void GPU::GPUSolver2::initializeDeviceData_byHostData() {
	// 将主机复制到设备
	if (!GlobalPara::basic::useGPU) {
		LogWriter::logAndPrintError("useGPU = false. @initializeDeviceData_byHostData.\n");
		exit(-1);
	}
	if (!hostMemoryAllocated) {
		LogWriter::logAndPrintError("Host data not allocated. @GPU::GPUSolver2::initializeDeviceData_byHostData()\n");
		exit(-1);
	}
	if (!deviceMemoryAllocated) {
		LogWriter::logAndPrintError("Device data not allocated. @GPU::GPUSolver2::initializeDeviceData_byHostData()\n");
		exit(-1);
	}
	if (deviceDataInitialized) {
		LogWriter::logAndPrintError("Device data has initialized. @GPU::GPUSolver2::initializeDeviceData_byHostData()\n");
		exit(-1);
	}

	host_to_device();// 复制host数据到device
	deviceDataInitialized = true;
}

void GPU::GPUSolver2::device_to_host() {
	// 将device端的数据复制到host端
	GPU::NodeSoA::cuda_memcpy(&node_host, &node_device, cudaMemcpyDeviceToHost);
	GPU::ElementSoA::cuda_memcpy(&element_host, &element_device, cudaMemcpyDeviceToHost);
	GPU::ElementFieldSoA::cuda_memcpy(&elementField_host, &elementField_device, cudaMemcpyDeviceToHost);
	GPU::EdgeSoA::cuda_memcpy(&edge_host, &edge_device, cudaMemcpyDeviceToHost);
	GPU::BoundarySetMap::cuda_memcpy(&boundary_host, &boundary_device, cudaMemcpyDeviceToHost);
}

void GPU::GPUSolver2::host_to_device() {
	// 将host端的数据复制到device端
	GPU::NodeSoA::cuda_memcpy(&node_device, &node_host, cudaMemcpyHostToDevice);
	GPU::ElementSoA::cuda_memcpy(&element_device, &element_host, cudaMemcpyHostToDevice);
	GPU::ElementFieldSoA::cuda_memcpy(&elementField_device, &elementField_host, cudaMemcpyHostToDevice);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);
	GPU::BoundarySetMap::cuda_memcpy(&boundary_device, &boundary_host, cudaMemcpyHostToDevice);
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
	FieldWriter::getInstance()->freeNodeFieldData(this->outputNodeField);
	this->edgeField_host.free();
	this->elementFieldVariable_dt_host.free();
	// 更新状态
	hostMemoryAllocated = false;
}

void GPU::GPUSolver2::freeDeviceMemory() {
	if (!deviceMemoryAllocated) {
		LogWriter::logAndPrintError("cannot free device memory.\n");
		exit(-1);
	}
	// 双向数据
	this->node_device.cuda_free();
	this->element_device.cuda_free();
	this->elementField_device.cuda_free();
	this->edge_device.cuda_free();
	this->boundary_device.cuda_free();

	// device 数据
	//delete this->infPara_device;
	this->edgeField_device.cuda_free();
	this->elementFieldVariable_dt_device.cuda_free();

	// 更新状态
	deviceMemoryAllocated = false;
}

void GPU::GPUSolver2::updateResidualVector() {
	U2NITS::Output::get_residual_host(elementField_host, residualVector, U2NITS::Output::normType_infinity);
}

inline void check_cudaError_t(cudaError_t error) {
	if (error != cudaSuccess) {
		LogWriter::logError(cudaGetErrorString(error));
		exit(-1);
	}
}

void GPU::GPUSolver2::setGPUDevice() {
	if (hasSetGPUDevice) {
		LogWriter::logAndPrintError("has set GPU device, cannot set twice. @GPU::GPUSolver2::setGPUDevice()\n");
		exit(-1);
	}

	int nDevice = 0;
	int iDevice = 0;
	cudaDeviceProp deviceProp;
	check_cudaError_t(cudaGetDeviceCount(&nDevice));
	check_cudaError_t(cudaSetDevice(iDevice));
	check_cudaError_t(cudaGetDeviceProperties(&deviceProp, iDevice));
	std::stringstream ss;
	ss << "device num: " << nDevice << "; ";
	ss << "active device: " << iDevice << ",";
	ss << "name: " << deviceProp.name << ",";
	ss << "compute capability: " << deviceProp.major << "." << deviceProp.minor << "\n";

	LogWriter::logAndPrint(ss.str());
	// 更新状态
	hasSetGPUDevice = true;
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
		myfloat dx = node_host.xy[0][tmpNodeID1] - node_host.xy[0][tmpNodeID0];
		myfloat dy = node_host.xy[1][tmpNodeID1] - node_host.xy[1][tmpNodeID0];
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

	//// boundarySet
	//for (int i = 0; i < boundaryManager.boundaries.size(); i++) {
	//	boundary_host.type[i] = boundaryManager.boundaries[i].type;
	//}
}


void GPU::GPUSolver2::updateOutputNodeField() {
	FieldWriter::getInstance()->update_nodeField();
}
