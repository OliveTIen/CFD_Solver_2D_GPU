#ifndef GPU_SOLVER_2_H
#define GPU_SOLVER_2_H

//#include "dataType/ElementDataPack.h"
//#include "dataType/ElementAdjacent.h"
#include "datatype/NodeSoA.h"
#include "dataType/ElementSoA.h"
#include "datatype/FieldSoA.h"
#include "dataType/EdgeSoA.h"

namespace GPU {
	class GPUSolver2 {
	public:
		GPU::NodeSoA node_host;
		GPU::NodeSoA node_device;
		GPU::ElementSoA element_host;
		GPU::ElementSoA element_device;
		GPU::EdgeSoA edge_host;
		GPU::EdgeSoA edge_device;

		GPU::FieldSoA elementField_host;
		GPU::FieldSoA elementField_device;
		REAL* element_vruvp[4]{};// host 由elementField_host的U转换得到
		REAL* element_U_old[4]{};// host
		GPU::OutputNodeFieldSoA outputNodeField;

	public:
		void initialze();
		// GPU iteration， 包含几个GPU kernel
		void iteration();

		void finalize();
		
		// 更新host数据 用U更新ruvp U_old
		void update_host_data(void* _pFVM2D_);
		// device数据更新到host
		void device_to_host();
	private:
		void initialize_nodeHost(void* _pFVM2D_, int num_node);
		// 初始化element和elementField的host数据
		void initialize_elementHost(void* _pFVM2D_, int num_element);

		void initialize_edgeHost(void* _pFVM2D_, int num_edge);
	};
}

#endif // GPU_SOLVER_2_H