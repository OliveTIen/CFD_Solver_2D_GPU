#ifndef BASIC_ALGORITHM_GPU_H
#define BASIC_ALGORITHM_GPU_H
#include "../gpu/dataType/Define.h"

namespace GPU {
	namespace Math {
		
#ifdef max 
#undef max // �����undef�Ļ����ᵼ���Զ���max��������ʧ��
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
	}
}

#endif // !BASIC_ALGORITHM_GPU_H
