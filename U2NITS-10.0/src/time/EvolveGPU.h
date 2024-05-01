#ifndef __EVOLVE_GPU_H__
#define __EVOLVE_GPU_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace GPU {
	namespace Time {
		void EvolveDevice(myfloat dt, int flag_timeAdvance, ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary, SDevicePara& para);
		// ·Ç¶¨³£
		void evolveSingleStep_device(myfloat dt, GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device,
			BoundarySetMap& boundary, SDevicePara& para);
	}
}



#endif