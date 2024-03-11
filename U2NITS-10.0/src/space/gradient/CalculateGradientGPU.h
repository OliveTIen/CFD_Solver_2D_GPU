#ifndef _CALCULATE_GRADIENT_2_CUH_
#define _CALCULATE_GRADIENT_2_CUH_

#include "../../gpu/datatype/DataType.h"

namespace GPU {
	void calculateGradient2(GPU::ElementSoA& element_device, GPU::FieldSoA elementField_device);

	__global__ void calculateGradientKernel2(GPU::ElementSoA& element_device, GPU::FieldSoA elementField_device);
}

#endif