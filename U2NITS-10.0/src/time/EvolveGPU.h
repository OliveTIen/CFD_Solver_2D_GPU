#ifndef __EVOLVE_GPU_H__
#define __EVOLVE_GPU_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace GPU {
	namespace Time {
		void EvolveDevice(REAL dt, int flag, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);
		void EvolveExplicitDevice(REAL dt, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);
	}
}

#endif