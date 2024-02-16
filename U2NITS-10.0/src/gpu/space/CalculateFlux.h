#ifndef _CALCULATE_FLUX_H_
#define _CALCULATE_FLUX_H_

#include "../dataType/ElementDataPack.h"
#include "../dataType/ElementAdjacent.h"
#include "../datatype/EdgeSoA.h"

namespace GPU {
	void calculateFlux(ElementDataPack& element_device,EdgeSoA& edge_device);
}

#endif