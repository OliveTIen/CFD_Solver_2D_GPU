#ifndef _CALCULATE_FLUX_2_H_
#define _CALCULATE_FLUX_2_H_

#include "../dataType/ElementSoA.h"
#include "../dataType/EdgeSoA.h"
#include "../dataType/FieldSoA.h"

namespace GPU {
	void calculateFlux2(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device);

	REAL calculateLambdaFlux(REAL edgeU[4], REAL edgeN[2], REAL gamma, REAL length);
}

#endif