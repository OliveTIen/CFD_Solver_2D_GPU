#ifndef _CALCULATE_GRADIENT_2_CUH_
#define _CALCULATE_GRADIENT_2_CUH_

#include "../../gpu/datatype/DataType.h"

namespace GPU {
	void calculateGradient_old(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device);

	__global__ void calculateGradientKernel_old(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device);

	namespace Space {
		namespace Gradient {
			void GradientLeastSquare(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device, GPU::NodeSoA& node_device);
			__global__ void GradientLeastSquareKernel(
				GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device, GPU::NodeSoA& node_device
			);
		}
	}

}

#endif