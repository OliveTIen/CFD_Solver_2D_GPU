#ifndef PHYSICAL_CONVERT_KERNEL_H
#define PHYSICAL_CONVERT_KERNEL_H

#include "../gpu/dataType/Define.h"

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
		inline void U2ruvwp_host(const REAL U[5], REAL ruvwp[5], REAL gamma) {
			//// [别删]原始版本，用于理解和推导，不要删了给自己添麻烦
			//typedef REAL real;
			//real rho = U[0];
			//real rhou = U[1];
			//real rhov = U[2];
			//real rhow = U[3];
			//real rhoE = U[4];// 能量方程中的守恒量
			//real u = rhou / rho;
			//real v = rhov / rho;
			//real w = rhow / rho;
			//real E = rhoE / rho;
			//real V2 = u * u + v * v + w * w;
			//real e = E - 0.5 * V2;
			//// e = p/rho/(gamma-1)
			//real p = e * rho * (gamma - 1);

			// 简化版本
			ruvwp[0] = U[0];
			ruvwp[1] = U[1] / U[0];
			ruvwp[2] = U[2] / U[0];
			ruvwp[3] = U[3] / U[0];
			REAL V2 = ruvwp[1] * ruvwp[1] + ruvwp[2] * ruvwp[2] + ruvwp[3] * ruvwp[3];
			ruvwp[4] = (U[4] / U[0] - 0.5 * V2) * U[0] * (gamma - 1);
		}
	}
}

namespace GPU {
	namespace Math {
		// 如果不undef的话，会导致自定义max函数编译失败
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

		__device__ inline real max(real a, real b) {
			return (a > b) ? a : b;
		}

		__device__ inline real min(real a, real b) {
			return (a < b) ? a : b;
		}

		__device__ inline real abs(real a) {
			return (a > 0) ? a : (-a);
		}

		__device__ inline void U2ruvp(const real U[4], real ruvp[4], real gamma) {
			ruvp[0] = U[0];
			ruvp[1] = U[1] / U[0];
			ruvp[2] = U[2] / U[0];
			ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
		}
		__device__ inline void ruvp2U(const real ruvp[4], real U[4], real gamma) {
			U[0] = ruvp[0];
			U[1] = ruvp[0] * ruvp[1];
			U[2] = ruvp[0] * ruvp[2];
			U[3] = ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);//rhoE
		}
	}
}

#endif