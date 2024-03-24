#ifndef __GPUMATHKERNEL_H__
#define __GPUMATHKERNEL_H__

#include "Common.h"
#include <iostream>

namespace GPU {
	namespace Math {
		namespace Matrix {
			//constexpr real EPSILON = U2NITS::Math::EPSILON;
			using U2NITS::Math::EPSILON;

			inline void printMatrix(long i, const long j, REAL* matrix) {
				for (long x = 0; x < i; x++) {
					for (long y = 0; y < j; y++) {
						std::cout << *(matrix + x * j + y) << " ";
					}
					std::cout << std::endl;
				}
			}

			// mat: ixj, res: jxi
			inline void __device__ transpose(int i, int j, REAL* mat, REAL* res) {
				for (int x = 0; x < i; x++) {
					for (int y = 0; y < j; y++) {
						*(res + y * i + x) = *(mat + x * j + y);
					}
				}
			}

			inline void __device__ inv_2x2(REAL(*mat)[2]) {
				REAL det, inv_det, inv_a, inv_b, inv_c, inv_d;
				det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
				if (det >= 0 && det < EPSILON) det = EPSILON;
				else if (det > -EPSILON)det = -EPSILON;
				inv_det = 1.0 / det;
				inv_a = mat[1][1] * inv_det;
				inv_b = -mat[0][1] * inv_det;
				inv_c = -mat[1][0] * inv_det;
				inv_d = mat[0][0] * inv_det;
				mat[0][0] = inv_a;
				mat[0][1] = inv_b;
				mat[1][0] = inv_c;
				mat[1][1] = inv_d;
			}

			inline void __device__ inv_2x2(REAL(*mat)[2], REAL(*res)[2]) {
				REAL det, inv_det;
				det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
				if (det >= 0 && det < EPSILON) det = EPSILON;
				else if (det > -EPSILON)det = -EPSILON;
				inv_det = 1.0 / det;
				res[0][0] = mat[1][1] * inv_det;
				res[0][1] = -mat[0][1] * inv_det;
				res[1][0] = -mat[1][0] * inv_det;
				res[1][1] = mat[0][0] * inv_det;
			}

			inline void __device__ inv_2x2(real* mat) {
				REAL det, inv_det, inv_a, inv_b, inv_c, inv_d;
				det = mat[0 * 2 + 0] * mat[1 * 2 + 1] - mat[0 * 2 + 1] * mat[1 * 2 + 0];
				if (det >= 0 && det < EPSILON) det = EPSILON;
				else if (det > -EPSILON)det = -EPSILON;
				inv_det = 1.0 / det;
				inv_a = mat[1 * 2 + 1] * inv_det;
				inv_b = -mat[0 * 2 + 1] * inv_det;
				inv_c = -mat[1 * 2 + 0] * inv_det;
				inv_d = mat[0 * 2 + 0] * inv_det;
				mat[0 * 2 + 0] = inv_a;
				mat[0 * 2 + 1] = inv_b;
				mat[1 * 2 + 0] = inv_c;
				mat[1 * 2 + 1] = inv_d;
			}

			// mat1: ixj, mat2: jxk, res: ixk
			inline void __device__ mul_ixj_jxk(int i, int j, int k, REAL* mat1, REAL* mat2, REAL* res) {
				for (int x = 0; x < i; x++) {
					for (int y = 0; y < k; y++) {
						*(res + x * k + y) = 0.0;
						for (int z = 0; z < j; z++) {
							// res[x][y] += mat1[x][z] * mat2[z][y];
							*(res + x * k + y) += *(mat1 + x * j + z) * *(mat2 + z * k + y);
						}
					}
				}
			}

			inline void __device__ div_matrix_by_scalar(int nRow, int nCol, real* mat, real scalar) {
				// ¾ØÕó³ýÒÔ±êÁ¿
				for (int i = 0; i < nRow; i++) {
					for (int j = 0; j < nCol; j++) {
						mat[i * nCol + j] /= scalar;
					}
				}
			}

		}
	}
}


#endif // __GPUMATHKERNEL_H__