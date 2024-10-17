#ifndef MATRIX_4XN_H
#define MATRIX_4XN_H
#include "DefineType.h"
#include "../Env.h"
namespace GPU {
	struct Matrix_4xn {
		const int nrow = 4;
		myint ncol = 0;
		myfloat* data[4]{};

		void cuda_alloc(myint _col) {
			ncol = _col;
			for (int i = 0; i < nrow; i++) {
				cudaMalloc(&data[i], ncol * sizeof(myfloat));
			}
		}

		void cuda_free() {
			for (int i = 0; i < nrow; i++) {
				cudaFree(data[i]);
			}
		}
	};

}

#endif // !MATRIX_4XN_H
