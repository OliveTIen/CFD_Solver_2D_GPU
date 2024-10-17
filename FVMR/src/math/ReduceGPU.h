#ifndef REDUCE_GPU_H
#define REDUCE_GPU_H
#include "../gpu/datatype/DefineType.h"

namespace GPU {
	namespace Math {
		enum ReduceType {
			reduceType_min,
			reduceType_max
		};


		// 规约运算，结果存入dev_output[0]。
		// n_dev_input: dev_input数组大小，必须为2的幂。dev_output数组大小至少为ceil(n/block_size)
		// p_func_device要求是有__device__标记的函数指针，例如__device__ func_bin_myfloat p_operator_min = operator_min; 目前如果直接传参有问题(cuda error 13)
		void reduce_device(const myint n_dev_input, myfloat* dev_input, myfloat* dev_output, bool debug_info, ReduceType reduceType);
	}
}

#endif // !REDUCE_GPU_H
