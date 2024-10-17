#include "FieldSoA.h"
#include "../Env.h"

void GPU::ElementFieldSoA::cuda_memcpy(ElementFieldSoA* dst, const ElementFieldSoA* src, cudaMemcpyKind kind) {
	myint _num = dst->num;
	for (int i = 0; i < 4; i++) {
		cudaMemcpy(dst->U[i], src->U[i], _num * sizeof(myfloat), kind);
		cudaMemcpy(dst->Ux[i], src->Ux[i], _num * sizeof(myfloat), kind);
		cudaMemcpy(dst->Uy[i], src->Uy[i], _num * sizeof(myfloat), kind);
		cudaMemcpy(dst->Flux[i], src->Flux[i], _num * sizeof(myfloat), kind);
		cudaMemcpy(dst->ruvp[i], src->ruvp[i], _num * sizeof(myfloat), kind);
		cudaMemcpy(dst->Uold[i], src->Uold[i], _num * sizeof(myfloat), kind);
	}
}


void GPU::EdgeFieldSoA::cuda_memcpy(EdgeFieldSoA* dst, const EdgeFieldSoA* src, cudaMemcpyKind kind) {
	myint num = dst->num_edge;
	for (int i = 0; i < 4; i++) {
		cudaMemcpy(dst->Flux[i], src->Flux[i], num * sizeof(myfloat), kind);
	}
}

void GPU::NodeFieldSoA::cuda_memcpy(NodeFieldSoA* dst, const NodeFieldSoA* src, cudaMemcpyKind kind) {
	myint num = dst->num_node;
	for (int i = 0; i < 4; i++) {
		cudaMemcpy(dst->ruvp[i], src->ruvp[i], num * sizeof(myfloat), kind);
	}
}

void GPU::ElementFieldVariable_dt::cuda_memcpy(ElementFieldVariable_dt* dst, const ElementFieldVariable_dt* src, cudaMemcpyKind kind) {
	myint num = dst->num_element;
	cudaMemcpy(dst->alphaC, src->alphaC, num * sizeof(myfloat), kind);
}

// 计算规约操作中dev_output数组长度

myint GPU::ReduceHelper::get_dev_output_length(myint n) {
	/*
	规约操作要求输入数组长度n补齐为2的幂次方，否则某一步规约会出现奇数，导致末尾元素没有参与规约
	对于求最小值，可以补较大的数
	*/
	myint block_threads = GPU::get_max_threads_per_block();// 每个block的线程数
	myint threads_needed = n / 2;// n个元素，第1次规约需要n/2个线程
	myint blocks = threads_needed / block_threads + (threads_needed % block_threads > 0 ? 1 : 0);
	return blocks;
}

// 获取大于等于original_num的最小的2的幂次方

myint GPU::ReduceHelper::get_next_power_2_number(myint n) {
	/*
	当输入2,147,483,647时，next_pow2最终会增大到1,073,741,824，然后变成-2,147,483,648，再乘以2变成0，导致死循环
	因此需要判断myint的最大值
	size of int: 4，因此最大=2^(4*8)-1=2,147,483,647
	size of long: 4
	size of int*: 8
	size of long long: 8
	*/
	myint next_pow2 = 2;
	while (next_pow2 < n) {
		next_pow2 *= 2;
		if (next_pow2 < 0) {
			printf("error: next pow2 is out of range," __FILE__);
			throw "out of range";
		}
	}
	return next_pow2;
}
