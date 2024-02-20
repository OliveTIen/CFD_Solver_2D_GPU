#ifndef _CALCULATE_FLUX_2_H_
#define _CALCULATE_FLUX_2_H_

#include "../dataType/ElementSoA.h"
#include "../datatype/EdgeSoA.h"

namespace GPU {
	void calculateFlux2(ElementSoA& element_device, EdgeSoA& edge_device);
}

#endif