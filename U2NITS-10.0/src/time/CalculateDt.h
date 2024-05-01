#ifndef _CALCULATE_DT_H_
#define _CALCULATE_DT_H_

#include "../math/Math.h"
#include "../solvers/GPUSolver2.h"
#include <iostream>

namespace U2NITS {
	namespace Time {
		myfloat calculateGlobalDt(
			myfloat t, myfloat T, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]
		);

		// ���㱾��ʱ�䲽���첽ģʽ��û�����ݾ�������Euler
		void calculateLocalTimeStep_async_Euler(
			myfloat& dt, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, myint iElement,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]
		);

	}
}

#endif