#ifndef _CALCULATE_DT_H_
#define _CALCULATE_DT_H_

#include "../math/Math.h"
#include "../solvers/GPUSolver2.h"
#include <iostream>

namespace U2NITS {
	namespace Time {
		myfloat get_global_dt_host(
			myfloat t, myfloat T, myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt elementFieldVariable_dt, 
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4], int physicsModel_equation
		);

		// 计算本地时间步，异步模式（没有数据竞争），Euler
		void calculateLocalTimeStep_async_Euler(
			myfloat& dt, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, myint iElement,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]
		);

	}
}

#endif