#ifndef PHYSICAL_KERNEL_GPU_H
#define PHYSICAL_KERNEL_GPU_H

#include "Common.h"

namespace GPU {
	namespace Math {

		__device__ inline void U2ruvp(const myfloat U[4], myfloat ruvp[4], myfloat gamma) {
			ruvp[0] = U[0];
			ruvp[1] = U[1] / U[0];
			ruvp[2] = U[2] / U[0];
			ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
		}
		__device__ inline void ruvp2U(const myfloat ruvp[4], myfloat U[4], myfloat gamma) {
			U[0] = ruvp[0];
			U[1] = ruvp[0] * ruvp[1];
			U[2] = ruvp[0] * ruvp[2];
			U[3] = ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);//rhoE
		}
	}
}


#endif