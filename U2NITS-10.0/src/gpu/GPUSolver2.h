#ifndef GPU_SOLVER_2_H
#define GPU_SOLVER_2_H

#include "dataType/ElementDataPack.h"
#include "dataType/ElementAdjacent.h"
#include "dataType/ElementSoA.h"
#include "dataType/EdgeSoA.h"

namespace GPU {
	class GPUSolver2 {
	public:
		GPU::ElementSoA element_host;
		GPU::ElementSoA element_device;
		GPU::EdgeSoA edge_host;
		GPU::EdgeSoA edge_device;

	public:
		void initialze();

		void iteration();

		void finalize();
		
	private:
		void initialize_elementHost(void* _pFVM2D_, int num_element);

		void initialize_edgeHost(void* _pFVM2D_, int num_edge);
	};
}

#endif // GPU_SOLVER_2_H