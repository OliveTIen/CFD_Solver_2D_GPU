#ifndef GPU_SOLVER_2_H
#define GPU_SOLVER_2_H

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace GPU {
	class GPUSolver2 {
	public:
		// host/device双向数据
		GPU::NodeSoA node_host;
		GPU::NodeSoA node_device;
		GPU::ElementSoA element_host;
		GPU::ElementSoA element_device;
		GPU::FieldSoA elementField_host;
		GPU::FieldSoA elementField_device;
		GPU::EdgeSoA edge_host;
		GPU::EdgeSoA edge_device;
		GPU::BoundarySetMap boundary_host;// 线性表，将boundary的setID映射到type
		GPU::BoundarySetMap boundary_device;

		// host数据
		REAL* element_vruvp[4]{};// 存储U转换为的ruvp，用于计算节点ruvp以及计算时间步长Dt
		REAL* element_U_old[4]{};// 存储上一步U用于计算残差
		std::map<int, int> edge_periodic_pair;// 存储周期边界的边界对
		GPU::OutputNodeFieldSoA outputNodeField;

		// device数据
		GPU::DGlobalPara* infPara_device;// 远场参数在GPU的副本，用于传参[RAII类型]


	public:
		// 申请资源、初始化GPUID、host复制到device
		void initialze();
		void iterationTotalCPU(REAL& t, REAL T);
		void iterationHalfGPU(REAL& t, REAL T);
		void iterationGPU(REAL& t, REAL T);
		// 更新节点场
		void updateOutputNodeField();
		// 更新残差
		void updateResidual();
		// 释放资源
		void finalize();
		
	private:
		void initialize_nodeHost(void* _pFVM2D_, int num_node);
		// 初始化element和elementField的host数据
		void initialize_elementHost(void* _pFVM2D_, int num_element);
		// 初始化edge的host数据
		void initialize_edgeHost(void* _pFVM2D_, int num_edge);
		// 初始化边界条件 将周期边界变为内部边界
		void initialize_boundary(void* _pFVM2D_);

		// 仅host端的更新 用U更新ruvp U_old
		void update_ruvp_Uold();
		// 计算时间步长
		REAL calculateDt(REAL t, REAL T);
		// GPU iteration， 包含几个GPU kernel
		void iterationDevice(REAL dt);
		// device数据更新到host
		void device_to_host();
		// 仅host iteration
		void iterationHost(REAL dt);
		// host数据更新到device
		void host_to_device();
	};
}

#endif // GPU_SOLVER_2_H