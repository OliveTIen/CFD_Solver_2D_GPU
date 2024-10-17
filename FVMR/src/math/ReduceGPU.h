#ifndef REDUCE_GPU_H
#define REDUCE_GPU_H
#include "../gpu/datatype/DefineType.h"

namespace GPU {
	namespace Math {
		enum ReduceType {
			reduceType_min,
			reduceType_max
		};


		// ��Լ���㣬�������dev_output[0]��
		// n_dev_input: dev_input�����С������Ϊ2���ݡ�dev_output�����С����Ϊceil(n/block_size)
		// p_func_deviceҪ������__device__��ǵĺ���ָ�룬����__device__ func_bin_myfloat p_operator_min = operator_min; Ŀǰ���ֱ�Ӵ���������(cuda error 13)
		void reduce_device(const myint n_dev_input, myfloat* dev_input, myfloat* dev_output, bool debug_info, ReduceType reduceType);
	}
}

#endif // !REDUCE_GPU_H
