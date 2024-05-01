#ifndef BASIC_ALGORITHM_H
#define BASIC_ALGORITHM_H
#include "../gpu/dataType/DefineType.h"
#include <cmath>
namespace U2NITS {
	namespace Math {
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
		inline myfloat max(myfloat a, myfloat b) {
			return (a > b) ? a : b;
		}

		inline myfloat min(myfloat a, myfloat b) {
			return (a < b) ? a : b;
		}

		inline myfloat abs(myfloat a) {
			return (a > 0) ? a : (-a);
		}

		inline myfloat distance2D(myfloat x1, myfloat y1, myfloat x2, myfloat y2) {
			myfloat dx = x2 - x1;
			myfloat dy = y2 - y1;
			return sqrt(dx * dx + dy * dy);
		}
		// 二维向量点积 u·v=u[0] * v[0] + u[1] * v[1]
		inline myfloat dot2D(myfloat* u, myfloat* v) {
			return u[0] * v[0] + u[1] * v[1];
		}
		// 向量数乘 v[i] *= scalar;
		inline void vector_times_scalar(myfloat* v, int v_length, myfloat scalar) {
			for (int i = 0; i < v_length; i++) {
				v[i] *= scalar;
			}
		}
		// 向量相加 result[i] = u[i] + v[i]
		inline void vector_add(const myfloat* u, const myfloat* v, myfloat* result, int length) {
			for (int i = 0; i < length; i++) {
				result[i] = u[i] + v[i];
			}
		}
		// 向量相加 v[i] += u[i];
		inline void vector_addto_v(const myfloat* u, myfloat* v, int length) {
			for (int i = 0; i < length; i++) {
				v[i] += u[i];
			}
		}
	}
}

#endif