#ifndef CALCULATE_DT_GPU_H
#define CALCULATE_DT_GPU_H
#include "../gpu/datatype/DataType.h"

namespace GPU {
	namespace Time {
		real calculateGlobalDt(
			real t, real T, real gamma, real Re, real Pr, real CFL, real Rcpcv,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, real* element_vruvp[4]
		);
	}
}

#endif // !CALCULATE_DT_GPU_H
