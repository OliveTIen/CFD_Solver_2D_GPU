#ifndef __EVOLVE_H__
#define __EVOLVE_H__

#include "../datatype/FieldSoA.h"
#include "../datatype/ElementSoA.h"
#include "../datatype/EdgeSoA.h"
#include <map>

namespace GPU {
	namespace Time {
		void EvolveHost(REAL dt, int flag, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);
		void EvolveExplicitHost(REAL dt, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);
	}
}

#endif