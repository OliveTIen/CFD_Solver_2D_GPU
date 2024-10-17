#include "ReduceGPU.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../math/BasicAlgorithmGPU.h"

/*
��2��kernel���ϲ�ȫ���ڴ�
ÿ���̴߳����������ݣ����������ݵĿռ���루������stride������ͬ�ģ����ƻ�������
*/
__global__ void reduce_device_kernel(myfloat* input, myfloat* output, unsigned int n, func_bin_myfloat p_func) {
    // Determine this thread's various ids
    unsigned int block_size = blockDim.x;
    unsigned int thread_id = threadIdx.x;
    unsigned int block_id = blockIdx.x;

    /*
    ÿ���̲߳���input[block_start]��input[block_start + stride]������input[block_start]
    ��˷�����ֱ�����ʣ��1����
    */
    unsigned int block_start = block_id * block_size * 2 + thread_id;// ��ǰ�̵߳��������
    for (unsigned int stride = block_size; stride > 0; stride /= 2) {
        if (thread_id < stride && // ��ǰ�̵߳���������ڷ�Χ�� 
            block_start + stride < n) // ��ǰ�̵߳��Ҳ������ڷ�Χ��
        {
            //input[block_start] += input[block_start + stride];

            p_func(input[block_start], input[block_start + stride]);
        }
        // ͬ�������߳�
        __syncthreads();
    }

    // �����ʣ�µ�1��������output
    if (!thread_id) {
        output[block_id] = input[block_start];
    }
}

inline void reduce_device_show_debug_info(unsigned int remaining, unsigned int blocks, unsigned int threads_needed) {
    printf("Launching kernels:\n");
    printf("remaining: %u\n", remaining);
    printf("blocks: %u\n", blocks);
    printf("threads_needed: %u\n", threads_needed);
    printf("\n");
}

void GPU::Math::reduce_device(const myint n, myfloat* dev_input, myfloat* dev_output, bool debug_info, ReduceType reduceType) {
    /*
    n dev_input���ܳ��ȣ�Ҫ��Ϊ2���ݡ�������ĳһ����nΪ���������һ��Ԫ��û�в������㣬���½����׼ȷ
    dev_input ����Լ���� ��СΪn
    dev_output �м����鼰����������� ��СΪ n/block_size ����ȡ���������Ѿ������
    */

    // ���n�Ƿ�Ϊ2���ݡ����� https://blog.csdn.net/qq_39360985/article/details/78628550
    if ((n & n - 1) == 0) {
        //printf("%d is pow of 2\n", n);
    }
    else {
        printf("warning: %d is NOT pow of 2\n", n);
    }

    // ���ݹ�Լ���ͣ�ѡ���Ӧ��˫Ŀ�����(����ָ��)
    func_bin_myfloat p_func_host;// ����ָ��
    switch (reduceType) {
    case reduceType_min:
        cudaMemcpyFromSymbol(&p_func_host, GPU::Math::p_operator_min, sizeof(func_bin_myfloat));// device to host
        break;
    case reduceType_max:
        cudaMemcpyFromSymbol(&p_func_host, GPU::Math::p_operator_max, sizeof(func_bin_myfloat));// device to host
        break;
    default:
        cudaMemcpyFromSymbol(&p_func_host, GPU::Math::p_operator_min, sizeof(func_bin_myfloat));
    }
    getLastCudaError("cudaMemcpyFromSymbol @ reduce_device failed.");

    const int block_threads = GPU::get_max_threads_per_block();// ÿ��block������߳���
    unsigned int threads_needed = n / 2; // �߳�����
    unsigned int blocks = threads_needed / block_threads + (threads_needed % block_threads > 0 ? 1 : 0); // block����
    unsigned int remaining = n; // ��Ҫ��͵�Ԫ������
    while (remaining > 1) {
        // ��ʾ������Ϣ
        if (debug_info) {
            reduce_device_show_debug_info(remaining, blocks, threads_needed);
        }

        // ���ú˺���
        reduce_device_kernel <<<blocks, block_threads>>> (dev_input, dev_output, remaining, p_func_host);

        // ������һ�������������Ϣ
        remaining = blocks; // ��һ����Ҫ��͵�Ԫ������
        threads_needed = remaining / 2; // ��һ����Ҫ���߳�����
        blocks = threads_needed / block_threads + (threads_needed % block_threads ? 1 : 0); // ��һ����Ҫ��block����

        // ����ָ�롣��û�д������ݵ�host�������޿���
        if (remaining > 1) {
            myfloat* dev_temp = dev_input;
            dev_input = dev_output;
            dev_output = dev_temp;
        }
    }
}
