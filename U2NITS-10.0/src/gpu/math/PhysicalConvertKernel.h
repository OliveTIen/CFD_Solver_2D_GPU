#ifndef PHYSICAL_CONVERT_KERNEL_H
#define PHYSICAL_CONVERT_KERNEL_H

#include "../dataType/Define.h"

namespace GPU {
	namespace Physics {
		// 默认标记为global，global是不是可以既可以在host也可以在device运行？
		// 将U数组批量转换为ruvp数组
		inline void U2ruvp_host(REAL* vU[4], REAL* vruvp[4], REAL gamma, long n) {

			//for (long i = 1; i < n; i++) {
			//	//U:	rho,rho_u,rho_v,rho_E
			//	//ruvp: rho,u,v,p
			//	REAL U[4]{ vU[0][i],vU[1][i],vU[2][i],vU[3][i] };
			//	REAL ruvp[4]{};

			//	ruvp[0] = U[0];
			//	ruvp[1] = U[1] / U[0];
			//	ruvp[2] = U[2] / U[0];
			//	ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);

			//	for (int j = 0; j < 4; j++) {
			//		vruvp[j][i] = ruvp[j];
			//	}
			//}
		}
		/*

		void Math_2D::U_2_ruvp(const double* U, double* ruvp, double gamma) {
			//U:	rho,rho_u,rho_v,rho_E
			//ruvp: rho,u,v,p
			ruvp[0] = U[0];
			ruvp[1] = U[1] / U[0];
			ruvp[2] = U[2] / U[0];
			ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
		}

		*/
	}
}


#endif