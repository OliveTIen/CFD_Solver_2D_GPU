#ifndef _CALCULATE_GRADIENT_2_H_
#define _CALCULATE_GRADIENT_2_H_

#include "../datatype/ElementSoA.h"
#include "../datatype/FieldSoA.h"

namespace GPU {
	void calculateGradient2(GPU::ElementSoA& element_device, GPU::FieldSoA elementField_device);

	__global__ void calculateGradientKernel2(GPU::ElementSoA& element_device, GPU::FieldSoA elementField_device);
}

#endif