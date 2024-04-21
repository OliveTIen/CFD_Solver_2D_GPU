
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
#include "../global/CExit.h"
#include "../time/CIteration.h"


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
	allocateMemory(num_node, num_element, num_edge, num_boundary);
}

void GPU::GPUSolver2::allocateMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
	allocateHostMemory(num_node, num_element, num_edge, num_boundary);
	if (GlobalPara::basic::useGPU) {
		setGPUDevice();
		allocateDeviceMemory(num_node, num_element, num_edge, num_boundary);
	}
}

void GPU::GPUSolver2::allocateHostMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
	if (hostMemoryAllocated) {
		LogWriter::logAndPrintError("Has allocated host memory.\n");
		CExit::saveAndExit(-1);
	}
	// ����host�ڴ�
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
	// ����״̬
	hostMemoryAllocated = true;
}

void GPU::GPUSolver2::allocateDeviceMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary) {
	if (deviceMemoryAllocated) {
		LogWriter::logAndPrintError("Has allocated device memory��cannot allocate twice. @GPUSolver2::allocateDeviceMemory\n");
		CExit::saveAndExit(-1);
	}


	if (!hasSetGPUDevice) {
		LogWriter::logAndPrintError("Hasn't set device. Cannot allocate device memory. @GPUSolver2::allocateDeviceMemory\n");
		exit(-1);
	}
	// ����device�ڴ�
	try {
		// ˫�򿽱�����
		node_device.cuda_alloc(num_node);
		element_device.cuda_alloc(num_element);
		elementField_device.cuda_alloc(num_element);
		edge_device.cuda_alloc(num_edge);
		boundary_device.cuda_alloc(num_boundary);

		// ������GPU������
		using namespace GlobalPara::boundaryCondition::_2D;// using�������򣬲��õ��Ļ���
		infPara_device = new GPU::DGlobalPara(inf::ruvp, inlet::ruvp, outlet::ruvp, 4
			, &GlobalPara::constant::R, &GlobalPara::constant::gamma, &GlobalPara::inviscid_flux_method::flux_conservation_scheme);
		sDevicePara.initialize(
			GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::c0, GlobalPara::constant::gamma, U2NITS::Math::EPSILON, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::constant::mu,
			GlobalPara::inviscid_flux_method::flag_reconstruct, GlobalPara::inviscid_flux_method::flag_gradient,
			GlobalPara::boundaryCondition::_2D::inf::ruvp, GlobalPara::boundaryCondition::_2D::inlet::ruvp, GlobalPara::boundaryCondition::_2D::outlet::ruvp,
			GlobalPara::inviscid_flux_method::flux_conservation_scheme, GlobalPara::inviscid_flux_method::flux_limiter
		);

		// ����״̬
		deviceMemoryAllocated = true;
	}
	catch (const char* e) {
		// ! �쳣����δ���
		std::cout << e << std::endl;
		cudaError_t error = cudaError_t::cudaErrorDeviceUninitialized;
		exit(-1);
	}
}

void GPU::GPUSolver2::initializeData_byOldData() {

	initializeHostData_byOldData();// �����ڴ棬��FVM_2D�����ݳ�ʼ��
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

	// ��FVM_2D�����ݳ�ʼ��host����
	U2NITS::OldDataConverter oldDataConverter(FVM_2D::getInstance(), this);
	oldDataConverter.Convert_FVM2D_to_HostData();

	hostDataInitialized = true;
}

void GPU::GPUSolver2::initializeDeviceData_byHostData() {
	// ���������Ƶ��豸
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

	host_to_device();// ����host���ݵ�device
	deviceDataInitialized = true;
}

void GPU::GPUSolver2::iteration(real& t, real T) {
	int useGPU = GlobalPara::basic::useGPU;
	if (useGPU == 0) {
		CIteration::iteration_host_2(t, T, this);
	}
	else if (useGPU == 1) {
		LogWriter::logAndPrintWarning("GPU mode is in development. Enter debugging mode. @GPU::GPUSolver2::iteration\n");
		CIteration::iteration_useGPU(t, T, this);
	}
	else {
		LogWriter::logAndPrintError("Invalid parameter: useGPU. @GPU::GPUSolver2::iteration\n");
		exit(-1);
	}

	if (!iterationStarted) {
		iterationStarted = true;
	}
}

void GPU::GPUSolver2::device_to_host() {
	// ��device�˵����ݸ��Ƶ�host��
	GPU::NodeSoA::cuda_memcpy(&node_host, &node_device, cudaMemcpyDeviceToHost);
	GPU::ElementSoA::cuda_memcpy(&element_host, &element_device, cudaMemcpyDeviceToHost);
	GPU::FieldSoA::cuda_memcpy(&elementField_host, &elementField_device, cudaMemcpyDeviceToHost);
	GPU::EdgeSoA::cuda_memcpy(&edge_host, &edge_device, cudaMemcpyDeviceToHost);

}

void GPU::GPUSolver2::host_to_device() {
	// ��host�˵����ݸ��Ƶ�device��
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
	// host host/device˫������
	this->node_host.free();
	this->element_host.free();
	this->elementField_host.free();
	this->edge_host.free();
	this->boundary_host.free();
	// host host����
	for (int j = 0; j < 4; j++) {
		delete[] element_vruvp[j];
		delete[] element_U_old[j];
	}
	this->outputNodeField.free();
	// ����״̬
	hostMemoryAllocated = false;
}

void GPU::GPUSolver2::freeDeviceMemory() {
	if (!deviceMemoryAllocated) {
		LogWriter::logAndPrintError("cannot free device memory.\n");
		exit(-1);
	}
	// device ˫������
	this->node_device.cuda_free();
	this->element_device.cuda_free();
	this->elementField_device.cuda_free();
	this->edge_device.cuda_free();
	this->boundary_device.cuda_free();

	// device device����
	delete this->infPara_device;
	// ����״̬
	deviceMemoryAllocated = false;
}

void GPU::GPUSolver2::updateResidual() {
	ResidualCalculator::cal_residual_GPU(element_U_old, elementField_host, ResidualCalculator::NORM_INF, residualVector);
}

void GPU::GPUSolver2::getResidual() {
	ResidualCalculator::get_residual_functionF(elementField_host, residualVector, ResidualCalculator::NORM_INF);
}

void GPU::GPUSolver2::setGPUDevice() {
	if (hasSetGPUDevice) {
		LogWriter::logAndPrintError("has set GPU device, cannot set twice. @GPU::GPUSolver2::setGPUDevice()\n");
		exit(-1);
	}

	int nDevice = 0;
	int iDevice = 0;	
	cudaError_t error;
	// ��ȡ�豸��
	error = cudaGetDeviceCount(&nDevice);
	if (error != cudaSuccess || nDevice <= 0) {
		LogWriter::log(cudaGetErrorString(error), LogWriter::Error);
		exit(-1);
	}
	// ���ù����豸
	error = cudaSetDevice(iDevice);
	if (error != cudaSuccess) {
		LogWriter::log(cudaGetErrorString(error), LogWriter::Error);
		exit(-1);
	}
	std::stringstream ss;
	ss << "Num of GPU: " << nDevice << ", " << "active GPU: " << iDevice << "\n";
	LogWriter::logAndPrint(ss.str());
	// ����״̬
	hasSetGPUDevice = true;
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


	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
	BoundaryManager_2D& boundaryManager = pFVM2D->boundaryManager;
	// ���ڱ߽�
	/**
	* ���ڱ߽��ж��ִ�����
	* һ������EdgeSoA������Աint periodicPair���ŵ��ǲ��ҿ죬ȱ����ռ�ÿռ��
	* ��һ������map��hashTable�洢���������ݽṹ��д�Ƚϸ��ӣ���GPU�����бȽ��鷳��
	* ��ΪҪ���±�д__device__�������ڲ���
	*
	* ĿǰGPU���õ�һ�ַ���
	*/
	for (VirtualBoundarySet_2D& boundary : boundaryManager.boundaries) {
		int bType = boundary.type;
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {
			// ��ĳ�߽缯set�����ڱ߽磬������edge
			for (Edge_2D* pEdge : boundary.pEdges) {
				int edgeID = pEdge->GPUID;
				int edgeID_pair = boundaryManager.get_pairEdge_periodic(pEdge)->GPUID;
				int elementL_pair = this->edge_host.elementL[edgeID_pair];
				// ��pair��elementL�����Լ���elementR
				this->edge_host.elementR[edgeID] = elementL_pair;
				// ��map�洢pair
				edge_periodic_pair.insert(std::pair<int, int>(edgeID, edgeID_pair));
				// ��EdgeSoA��Ա�洢pair
				this->edge_host.periodicPair[edgeID] = edgeID_pair;
			}
		}
	}

	// boundarySet
	for (int i = 0; i < boundaryManager.boundaries.size(); i++) {
		boundary_host.type[i] = boundaryManager.boundaries[i].type;
	}
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
