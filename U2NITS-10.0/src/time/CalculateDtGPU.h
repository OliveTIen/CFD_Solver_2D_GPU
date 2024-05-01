#ifndef CALCULATE_DT_GPU_H
#define CALCULATE_DT_GPU_H
#include "../gpu/datatype/DataType.h"

namespace GPU {
	namespace Time {
		myfloat calculateGlobalDt(
			myfloat t, myfloat T, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]
		);
	}
}

#endif // !CALCULATE_DT_GPU_H
