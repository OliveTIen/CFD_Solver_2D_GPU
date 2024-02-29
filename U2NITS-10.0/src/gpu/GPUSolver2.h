#ifndef GPU_SOLVER_2_H
#define GPU_SOLVER_2_H

//#include "dataType/ElementDataPack.h"
//#include "dataType/ElementAdjacent.h"
#include "datatype/NodeSoA.h"
#include "dataType/ElementSoA.h"
#include "datatype/FieldSoA.h"
#include "dataType/EdgeSoA.h"
#include <map>

namespace GPU {
	class GPUSolver2 {
	public:
		GPU::NodeSoA node_host;
		GPU::NodeSoA node_device;
		GPU::ElementSoA element_host;
		GPU::ElementSoA element_device;
		GPU::FieldSoA elementField_host;
		GPU::FieldSoA elementField_device;
		GPU::EdgeSoA edge_host;
		GPU::EdgeSoA edge_device;

		REAL* element_vruvp[4]{};// 存储U转换为的ruvp，用于计算节点ruvp以及计算时间步长Dt
		REAL* element_U_old[4]{};// 存储上一步U用于计算残差
		std::map<int, int> edge_periodic_pair;// 存储周期边界的边界对
		GPU::OutputNodeFieldSoA outputNodeField;
		//GPU::PeriodicBoundary periodicBoundary;// 仅临时用途，以后会优化掉 
		/* 占据大量空间，虽然查找复杂度为O(1)，但有许多0元素
		*/

	public:
		// 申请资源、初始化GPUID、host复制到device
		void initialze();
		void iteration(REAL& t, REAL T);
		// 更新节点场
		void updateOutputNodeField();
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