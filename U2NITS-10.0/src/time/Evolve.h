#ifndef __EVOLVE_H__
#define __EVOLVE_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace U2NITS {
	namespace Time {
		using namespace GPU;
		void EvolveHost(REAL dt, int flag, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);
		void EvolveExplicitHost(REAL dt, ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);
	}
}


#endif