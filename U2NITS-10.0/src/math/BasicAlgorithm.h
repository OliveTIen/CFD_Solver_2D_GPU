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
		inline real max(real a, real b) {
			return (a > b) ? a : b;
		}

		inline real min(real a, real b) {
			return (a < b) ? a : b;
		}

		inline real abs(real a) {
			return (a > 0) ? a : (-a);
		}

		inline real distance2D(real x1, real y1, real x2, real y2) {
			real dx = x2 - x1;
			real dy = y2 - y1;
			return sqrt(dx * dx + dy * dy);
		}
		// 二维向量点积 u·v=u[0] * v[0] + u[1] * v[1]
		inline real dot2D(real* u, real* v) {
			return u[0] * v[0] + u[1] * v[1];
		}
		// 向量数乘 v[i] *= scalar;
		inline void vector_times_scalar(real* v, int v_length, real scalar) {
			for (int i = 0; i < v_length; i++) {
				v[i] *= scalar;
			}
		}
		// 向量相加 result[i] = u[i] + v[i]
		inline void vector_add(const real* u, const real* v, real* result, int length) {
			for (int i = 0; i < length; i++) {
				result[i] = u[i] + v[i];
			}
		}
		// 向量相加 v[i] += u[i];
		inline void vector_addto_v(const real* u, real* v, int length) {
			for (int i = 0; i < length; i++) {
				v[i] += u[i];
			}
		}
	}
}

#endif