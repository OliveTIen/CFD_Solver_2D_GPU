#include "ReduceGPU.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../math/BasicAlgorithmGPU.h"

/*
第2版kernel，合并全局内存
每个线程处理两个数据，这两个数据的空间距离（即步长stride）是相同的，类似滑动窗口
*/
__global__ void reduce_device_kernel_1(float* input, float* output, unsigned int n, func_bin_myfloat p_func) {
    // Determine this thread's various ids
    unsigned int block_size = blockDim.x;
    unsigned int thread_id = threadIdx.x;
    unsigned int block_id = blockIdx.x;

    /*
    每个线程操作input[block_start]和input[block_start + stride]，存入input[block_start]
    如此反复，直到最后剩下1个数
    */
    unsigned int block_start = block_id * block_size * 2 + thread_id;// 当前线程的左操作数
    for (unsigned int stride = block_size; stride > 0; stride /= 2) {
        if (thread_id < stride && // 当前线程的左操作数在范围内 
            block_start + stride < n) // 当前线程的右操作数在范围内
        {
            //input[block_start] += input[block_start + stride];

            p_func(input[block_start], input[block_start + stride]);
        }
        // 同步所有线程
        __syncthreads();
    }

    // 将最后剩下的1个数存入output
    if (!thread_id) {
        output[block_id] = input[block_start];
    }
}

//__device__ func_bin_myfloat p_operator_min = GPU::Math::operator_min;

void GPU::Math::reduce_device(const myint n, myfloat* dev_input, myfloat* dev_output, bool debug_info, ReduceType reduceType) {
    /*
   n dev_input的总长度，要求为2的幂。否则在某一步，当n为奇数时，最后一个元素没有参与运算
   dev_input 待规约数组 大小为n
   dev_output 中间数组及最终输出数组 大小为 n/block_size 向上取整。事先已经分配好
   */

    // 检查n是否为2的幂。参照 https://blog.csdn.net/qq_39360985/article/details/78628550
    if ((n & n - 1) == 0) {
        //printf("%d is pow of 2\n", n);
    }
    else {
        printf("warning: %d is NOT pow of 2\n", n);
    }

    func_bin_myfloat p_func_host;// reduce kernel所用到的函数的指针

    switch (reduceType) {
    case reduceType_min:
        cudaMemcpyFromSymbol(&p_func_host, GPU::Math::p_operator_min, sizeof(func_bin_myfloat));// devic to host
        break;
    default:
        cudaMemcpyFromSymbol(&p_func_host, GPU::Math::p_operator_min, sizeof(func_bin_myfloat));
    }

    getLastCudaError("cudaMemcpyFromSymbol, reduce_device failed.");

    const int block_threads = GPU::get_max_threads_per_block();
    unsigned int threads_needed = n / 2; // we'll need one thread to add every 2 elements
    unsigned int blocks = threads_needed / block_threads +  // we'll need this many blocks
        (threads_needed % block_threads > 0 ? 1 : 0); // plus one extra if threads_needed

    unsigned int remaining = n; // tracks number of elements left to add
    while (remaining > 1) {
        if (debug_info) {
            printf("Launching kernels:\n");
            printf("remaining: %u\n", remaining);
            printf("blocks: %u\n", blocks);
            printf("threads_needed: %u\n", threads_needed);
            printf("\n");
        }

        // call the kernel
        reduce_device_kernel_1 <<<blocks, block_threads>>> (dev_input, dev_output, remaining, p_func_host);

        // re-compute our size information for the next iteration
        remaining = blocks; // After the previous kernel call, each block has reduced its chunk down to a single partial sum
        threads_needed = remaining / 2; // each thread added 2 elements
        blocks = threads_needed / block_threads + (threads_needed % block_threads ? 1 : 0); // again, might need one extra block if threads_needed
        // is not evenly divisible by block_threads

        // 交换指针。并没有传递数据到host，几乎无开销
        if (remaining > 1) {
            float* dev_temp = dev_input;
            dev_input = dev_output;
            dev_output = dev_temp;
        }
    }
}
