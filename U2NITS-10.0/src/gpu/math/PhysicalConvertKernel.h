#ifndef PHYSICAL_CONVERT_KERNEL_H
#define PHYSICAL_CONVERT_KERNEL_H

#include "../dataType/Define.h"

namespace GPU {
	namespace Physics {
		// 默认标记为global，global是不是可以既可以在host也可以在device运行？
		// 将U数组转换为ruvp数组
		inline void U2ruvp_host(const REAL U[4], REAL ruvp[4], REAL gamma) {
			// 参照 void Math_2D::U_2_ruvp(const double* U, double* ruvp, double gamma)
			// U:	rho,rho_u,rho_v,rho_E
			// ruvp: rho,u,v,p
			ruvp[0] = U[0];
			ruvp[1] = U[1] / U[0];
			ruvp[2] = U[2] / U[0];
			ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
		}

		inline void ruvp2U_host(const REAL ruvp[4], REAL U[4], REAL gamma) {
			// 参照 void Math_2D::ruvp_2_U(const double* ruvp, double* U, double gamma)
			// ruvp: rho,u,v,p
			// U:	rho,rho_u,rho_v,rho_E
			U[0] = ruvp[0];
			U[1] = ruvp[0] * ruvp[1];
			U[2] = ruvp[0] * ruvp[2];
			U[3] = ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);//rhoE
		}
	}
}


#endif