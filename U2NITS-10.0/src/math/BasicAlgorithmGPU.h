#ifndef BASIC_ALGORITHM_GPU_H
#define BASIC_ALGORITHM_GPU_H
#include "../gpu/dataType/DefineType.h"
#include "../gpu/Env.h"

namespace GPU {
	namespace Math {
		
#ifdef max 
#undef max // 如果不undef的话，会导致自定义max函数编译失败
#endif
#ifdef min
#undef min
#endif

		__device__ inline myfloat max(myfloat a, myfloat b) {
			return (a > b) ? a : b;
		}

		__device__ inline myfloat min(myfloat a, myfloat b) {
			return (a < b) ? a : b;
		}

		__device__ inline myfloat abs(myfloat a) {
			return (a > 0) ? a : (-a);
		}

		__device__ inline myfloat distance2D(myfloat x1, myfloat y1, myfloat x2, myfloat y2) {
			myfloat dx = x2 - x1;
			myfloat dy = y2 - y1;
			return sqrt(dx * dx + dy * dy);
		}
	}
}

#endif // !BASIC_ALGORITHM_GPU_H
