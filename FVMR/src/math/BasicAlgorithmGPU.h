#ifndef BASIC_ALGORITHM_GPU_H
#define BASIC_ALGORITHM_GPU_H
#include "../gpu/dataType/DefineType.h"
#include "../gpu/Env.h"

namespace GPU {
	namespace Math {
		/*
		inline functions
		*/

		
#ifdef max 
#undef max // 如果不undef的话，会导致自定义max函数编译失败
#endif
#ifdef min
#undef min
#endif

		__host__ __device__ inline myfloat max(myfloat a, myfloat b) {
			return (a > b) ? a : b;
		}

		__host__ __device__ inline myfloat min(myfloat a, myfloat b) {
			return (a < b) ? a : b;
		}

		__host__ __device__ inline myfloat abs(myfloat a) {
			return (a > 0) ? a : (-a);
		}

		__host__ __device__ inline myfloat distance2D(myfloat x1, myfloat y1, myfloat x2, myfloat y2) {
			myfloat dx = x2 - x1;
			myfloat dy = y2 - y1;
			return sqrt(dx * dx + dy * dy);
		}
		// 结果存入a
		__host__ __device__ inline void operator_min(myfloat& a, const myfloat& b) {
			if (a > b)a = b;
		}
		__host__ __device__ inline void operator_max(myfloat& a, const myfloat& b) {
			if (a < b)a = b;
		}
		// device函数指针
		extern __device__ func_bin_myfloat p_operator_min;
		extern __device__ func_bin_myfloat p_operator_max;

		/*
		kernels
		*/
		// 带权的向量除法 v1 /= v2 .* weight
		__global__ void vector_weighted_divide_kernel(myint length, myfloat* v1, const myfloat* v2, myfloat weight);
		// 带权的向量倒数 v1 = (v2 .* weight) / v1
		__global__ void vector_weighted_reciprocal_kernel(myint length, myfloat* v1, const myfloat* v2, myfloat weight);
		// 带权的向量加法 v1 += v2 .* weight
		__global__ void vector_weighted_add_kernel(myint length, myfloat* v1, const myfloat* v2, myfloat weight);
		// 向量点积 v1 += v2 .* v3
		__global__ void vector_dot_product_add_kernel(myint v_size, myfloat* v1, const myfloat* v2, const myfloat* v3);

		/*
		functions, vector
		*/
		// 给数组指定范围[start, end)的元素赋值
		void assign_elements_in_array_device(myint start, myint end, myfloat* arr_dev, myfloat value);
	}
}

#endif // !BASIC_ALGORITHM_GPU_H
