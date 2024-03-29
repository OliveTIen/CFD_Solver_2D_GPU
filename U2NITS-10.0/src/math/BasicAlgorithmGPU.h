#ifndef BASIC_ALGORITHM_GPU_H
#define BASIC_ALGORITHM_GPU_H
#include "../gpu/dataType/Define.h"

namespace GPU {
	namespace Math {
		
#ifdef max 
#undef max // 如果不undef的话，会导致自定义max函数编译失败
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

		__device__ inline real distance2D(real x1, real y1, real x2, real y2) {
			real dx = x2 - x1;
			real dy = y2 - y1;
			return sqrt(dx * dx + dy * dy);
		}
	}
}

#endif // !BASIC_ALGORITHM_GPU_H
