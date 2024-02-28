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
		REAL* element_vruvp[4]{};// host ��elementField_host��Uת���õ�
		REAL* element_U_old[4]{};// host
		GPU::OutputNodeFieldSoA outputNodeField;

	public:
		void initialze();
		// GPU iteration�� ��������GPU kernel
		void iteration();

		void finalize();
		
		// ����host���� ��U����ruvp U_old
		void update_host_data(void* _pFVM2D_);
		// device���ݸ��µ�host
		void device_to_host();
	private:
		void initialize_nodeHost(void* _pFVM2D_, int num_node);
		// ��ʼ��element��elementField��host����
		void initialize_elementHost(void* _pFVM2D_, int num_element);

		void initialize_edgeHost(void* _pFVM2D_, int num_edge);
	};
}

#endif // GPU_SOLVER_2_H