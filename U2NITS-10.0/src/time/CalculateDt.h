#ifndef _CALCULATE_DT_H_
#define _CALCULATE_DT_H_

#include "../math/Math.h"
#include "../solvers/GPUSolver2.h"
#include <iostream>

namespace U2NITS {
	namespace Time {
		myfloat get_global_dt_host(
			myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt, 
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation
		);

		myfloat get_global_dt_host_0529beta(
			myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation
		);

		// ���㱾��ʱ�䲽���첽ģʽ��û�����ݾ�������Euler
		void calculateLocalTimeStep_async_Euler_kernel(
			myfloat& dt, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, myint iElement,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]
		);

		// ����ֲ�ʱ�䲽��������λ��elementFieldVariable_dt.alphaC����
		void get_local_dt_host(
			myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt,
			GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation
		);

	}
}

#endif