#ifndef PHYSICAL_KERNEL_H
#define PHYSICAL_KERNEL_H

#include "CommonValue.h"
#include "StandardHeader.h"

namespace U2NITS {
	namespace Math {

		// 将U数组转换为ruvp数组
		inline void U2ruvp_host(const myfloat U[4], myfloat ruvp[4], myfloat gamma) {
			// U:	rho,rho_u,rho_v,rho_E 注意E=e+0.5*V^2
			// ruvp: rho,u,v,p
			ruvp[0] = U[0];
			ruvp[1] = U[1] / U[0];
			ruvp[2] = U[2] / U[0];
			ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
		}

		inline void ruvp2U_host(const myfloat ruvp[4], myfloat U[4], myfloat gamma) {
			// ruvp: rho,u,v,p
			// U:	rho,rho_u,rho_v,rho_E
			U[0] = ruvp[0];
			U[1] = ruvp[0] * ruvp[1];
			U[2] = ruvp[0] * ruvp[2];
			U[3] = ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);//rhoE
		}

		// 三维 守恒量转场变量
		inline void U2ruvwp_host_3d(const myfloat U[5], myfloat ruvwp[5], myfloat gamma) {

			/*
			e=Cv*t, Cv:R:Cp=1:(gamma-1):gama => R*t = (gamma-1)*Cv*t = (gamma-1)*e
			p=rho*R*t => p=rho*(gamma-1)*e = rho*(gamma-1)*(E-0.5*V2)
			*/
			ruvwp[0] = U[0];// rho
			ruvwp[1] = U[1] / U[0];// u
			ruvwp[2] = U[2] / U[0];// v
			ruvwp[3] = U[3] / U[0];// w
			myfloat V2 = ruvwp[1] * ruvwp[1] + ruvwp[2] * ruvwp[2] + ruvwp[3] * ruvwp[3];// 速度平方
			ruvwp[4] = (U[4] / U[0] - 0.5 * V2) * U[0] * (gamma - 1);// p
		}

		// 某方向上的马赫数
		inline myfloat getMach(const myfloat U[4], myfloat nx, myfloat ny, myfloat gamma) {
			myfloat ruvp[4]{};
			U2ruvp_host(U, ruvp, gamma);
			myfloat rho = ruvp[0];
			myfloat u = ruvp[1];
			myfloat v = ruvp[2];
			myfloat p = ruvp[3];
			myfloat un = u * nx + v * ny;
			myfloat a2 = gamma * p / rho;
			myfloat Ma2 = (u * u + v * v) / a2;
			return sqrt(Ma2);
		}

		
		// 批量转换
		void U_to_ruvp_in_batch(const myfloat* U[4], myfloat* ruvp[4], int length, myfloat gamma);
	}
}


#endif