#ifndef __EVOLVE_GPU_H__
#define __EVOLVE_GPU_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace GPU {
	namespace Time {
		void EvolveDevice(REAL dt, int flag_timeAdvance, ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary, SDevicePara& para);
		void EvolveExplicitDevice(REAL dt, ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary, SDevicePara& para);
		
		void TimeIntegration(REAL dt, ElementSoA& element_device, FieldSoA& elementField_device);
		__global__ void TimeIntegrationKernel(DReal& dt, ElementSoA& element, FieldSoA& elementField);
	}
}

#endif