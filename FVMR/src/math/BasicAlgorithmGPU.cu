#include "BasicAlgorithmGPU.h"
#include "../gpu/GPUGlobalFunction.h"

__global__ void GPU::Math::vector_weighted_divide_kernel(myint length, myfloat* v1, const myfloat* v2, myfloat weight) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < length) {
        v1[i] /= weight * v2[i];// �������ȼ� �����Ҳ�˷�������*=��/=
    }
}

__global__ void GPU::Math::vector_weighted_reciprocal_kernel(myint length, myfloat* v1, const myfloat* v2, myfloat weight) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < length) {
        v1[i] = weight * v2[i] / v1[i];
    }
}

__global__  void GPU::Math::vector_weighted_add_kernel(myint length, myfloat* v1, const myfloat* v2, myfloat weight) {
    myint i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < length) {
        v1[i] += weight * v2[i];
    }
}

__global__  void GPU::Math::vector_dot_product_add_kernel(myint v_size, myfloat* v1, const myfloat* v2, const myfloat* v3) {
    myint i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < v_size) {
        v1[i] += v2[i] * v3[i];
    }
}

__global__ void assign_elements_in_array_device_kernel(myint start, myint end, myfloat* arr, myfloat value) {
    const myint num = end - start;
    const myint idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        arr[start + idx] = value;
    }
}

void GPU::Math::assign_elements_in_array_device(myint start, myint end, myfloat* arr_dev, myfloat value) {
    // ������ָ����Χ[start, end)��Ԫ�ظ�ֵ
    myint num = end - start;
    if (num == 0)return;
    int block_size = GPU::get_max_threads_per_block();
	int grid_size = (num + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	assign_elements_in_array_device_kernel <<<grid, block>>> (start, end, arr_dev, value);
}

// device����ָ��
__device__ func_bin_myfloat GPU::Math::p_operator_min = operator_min;
__device__ func_bin_myfloat GPU::Math::p_operator_max = operator_max;