#ifndef _CALCULATE_DT_H_
#define _CALCULATE_DT_H_

#include "../math/Math.h"
#include "../solvers/GPUSolver2.h"
#include <iostream>

namespace U2NITS {

	namespace Time {
		real calculateGlobalDt(
			real t, real T, real gamma, real Re, real Pr, real CFL, real Rcpcv, 
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, real* element_vruvp[4]
		);

		// 计算本地时间步，异步模式（没有数据竞争），Euler
		void calculateLocalTimeStep_async_Euler(
			real& dt, real gamma, real Re, real Pr, real CFL, real Rcpcv, integer iElement,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, real* element_vruvp[4]
		);

	}
}

#endif