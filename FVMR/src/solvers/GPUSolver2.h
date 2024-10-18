#ifndef GPU_SOLVER_2_H
#define GPU_SOLVER_2_H

#include "../gpu/datatype/Datatype.h"
#include <map>
#include "../gpu/datatype/BoundaryV2.h"

namespace GPU {
	class GPUSolver2 {
	private:
		bool hostMemoryAllocated = false;
		bool hasSetGPUDevice = false;// 设备存在且已就绪
		bool deviceMemoryAllocated = false;
		bool iterationStarted = false;
		bool hostDataInitialized = false;
		bool deviceDataInitialized = false;
		static GPUSolver2* pInstance;// 实例的指针

	public:
		/*
		数据分为以下几种类型
		host文件输入数据。host和device情况都会使用
		host文件输出数据。host和device情况都会使用
		host运算数据。仅在host情况使用
		device运算数据。仅在device情况使用

		element_vruvp、element_U_old

		*/

		// host/device成对数据 若互相关联，一般用cuda_memcpy在初始化时从host复制到device
		GPU::NodeSoA node_host;
		GPU::NodeSoA node_device;
		GPU::ElementSoA element_host;
		GPU::ElementSoA element_device;
		GPU::ElementFieldSoA elementField_host;
		GPU::ElementFieldSoA elementField_device;
		GPU::ElementFieldVariable_dt elementFieldVariable_dt_host;// 互不关联，无需cuda_memcpy。仅作为中间变量
		GPU::ElementFieldVariable_dt elementFieldVariable_dt_device;
		GPU::EdgeSoA edge_host;
		GPU::EdgeSoA edge_device;
		GPU::EdgeFieldSoA edgeField_host;// 没有添加进host_to_device和device_to_host
		GPU::EdgeFieldSoA edgeField_device;// 没有添加进host_to_device和device_to_host
		GPU::BoundarySetMap boundary_host;// 旧版，存储了从set ID到set type的映射
		GPU::BoundarySetMap boundary_device;
		GPU::BoundaryV2 boundary_host_new;// 新版，存储了从set type到edge IDs的映射
		
		// host数据[未使用]
		myfloat* element_vruvp[4]{};// 存储U转换为的ruvp，用于计算节点ruvp以及计算时间步长Dt
		myfloat* element_U_old[4]{};// 存储上一步U用于计算残差, 4xn矩阵

		// host数据
		const int residualVectorSize = 4;
		myfloat residualVector[4]{ 1,1,1,1 };// 残差，residualVectorSize x 1向量
		std::map<int, int> edge_periodic_pair;// 存储周期边界的边界对
		GPU::OutputNodeFieldSoA outputNodeField;// nodeField_host

		// device数据[未使用]
		//GPU::DGlobalPara* infPara_device;// 已弃用
		SDevicePara sDevicePara;

	public:
		GPUSolver2();
		~GPUSolver2();
		// 为减少耦合，请使用SolverDataGetter间接调用该方法
		static GPUSolver2* getInstance();
		// should be called before allocating device memory
		void setGPUDevice();
		// 申请内存。应在读取field file后使用，因为要确定有多少单元、需要多大内存
		void allocateMemory();
		// including host and device
		void initializeData();
		// 更新节点场
		void updateOutputNodeField();
		// 获取残差
		void updateResidualVector();
		// 释放资源
		void freeMemory();

		void setIterationStarted(bool b) { iterationStarted = b; }
		bool isIterationStarted() { return iterationStarted; }
		bool isHostMemoryAllocated() { return hostMemoryAllocated; }
		bool isDeviceMemoryAllocated() { return deviceMemoryAllocated; }// first allocate memory
		bool isDeviceDataInitialized() { return deviceDataInitialized; }// then initialize data

		// device数据更新到host
		void device_to_host();
		// host数据更新到device
		void host_to_device();
	private:
		// using legacy FVM_2D data
		void initializeHostData_byOldData();
		void initializeOtherHostData();
		void initializeDeviceData_byHostData();
		void initialize_nodeHost(void* _pFVM2D_, int num_node);
		// 初始化element和elementField的host数据
		void initialize_elementHost(void* _pFVM2D_, int num_element);
		// 初始化edge的host数据
		void initialize_edgeHost(void* _pFVM2D_, int num_edge);
		// 初始化边界条件 将周期边界变为内部边界
		void initialize_boundary(void* _pFVM2D_);


		//void allocateMemory_kernel(const int num_node, const int num_element, const int num_edge, const int num_boundary);
		void allocateHostMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary);
		void allocateDeviceMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary);// 需要在setGPUDevice之后
		void freeHostMemory();
		void freeDeviceMemory();

	};
}

#endif // GPU_SOLVER_2_H