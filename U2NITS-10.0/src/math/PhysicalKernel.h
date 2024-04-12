#ifndef PHYSICAL_KERNEL_H
#define PHYSICAL_KERNEL_H

#include "Common.h"

namespace U2NITS {
	namespace Math {

		// 默认标记为global，global是不是可以既可以在host也可以在device运行？
		// 将U数组转换为ruvp数组
		inline void U2ruvp_host(const REAL U[4], REAL ruvp[4], REAL gamma) {
			// U:	rho,rho_u,rho_v,rho_E 注意E=e+0.5*V^2
			// ruvp: rho,u,v,p
			ruvp[0] = U[0];
			ruvp[1] = U[1] / U[0];
			ruvp[2] = U[2] / U[0];
			ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
		}

		inline void ruvp2U_host(const REAL ruvp[4], REAL U[4], REAL gamma) {
			// ruvp: rho,u,v,p
			// U:	rho,rho_u,rho_v,rho_E
			U[0] = ruvp[0];
			U[1] = ruvp[0] * ruvp[1];
			U[2] = ruvp[0] * ruvp[2];
			U[3] = ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);//rhoE
		}

		// 三维 守恒量转场变量
		inline void U2ruvwp_host_3d(const REAL U[5], REAL ruvwp[5], REAL gamma) {

			/*
			// [别删]原始版本，用于理解和推导，不要删了给自己添麻烦
			typedef REAL real;
			real rho = U[0];
			real rhou = U[1];
			real rhov = U[2];
			real rhow = U[3];
			real rhoE = U[4];// 能量方程中的守恒量
			real u = rhou / rho;
			real v = rhov / rho;
			real w = rhow / rho;
			real E = rhoE / rho;
			real V2 = u * u + v * v + w * w;
			real e = E - 0.5 * V2;
			// e = p/rho/(gamma-1)
			real p = e * rho * (gamma - 1);
			*/

			// 简化版本
			ruvwp[0] = U[0];
			ruvwp[1] = U[1] / U[0];
			ruvwp[2] = U[2] / U[0];
			ruvwp[3] = U[3] / U[0];
			REAL V2 = ruvwp[1] * ruvwp[1] + ruvwp[2] * ruvwp[2] + ruvwp[3] * ruvwp[3];
			ruvwp[4] = (U[4] / U[0] - 0.5 * V2) * U[0] * (gamma - 1);
		}

		// 某方向上的马赫数
		inline real getMach(const real U[4], real nx, real ny, real gamma) {
			real ruvp[4]{};
			U2ruvp_host(U, ruvp, gamma);
			real rho = ruvp[0];
			real u = ruvp[1];
			real v = ruvp[2];
			real p = ruvp[3];
			real un = u * nx + v * ny;
			real a2 = gamma * p / rho;
			real Ma2 = (u * u + v * v) / a2;
			return sqrt(Ma2);
		}
	}
}


#endif