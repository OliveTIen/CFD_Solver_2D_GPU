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

// �����Լ������dev_output���鳤��

myint GPU::ReduceHelper::get_dev_output_length(myint n) {
	/*
	��Լ����Ҫ���������鳤��n����Ϊ2���ݴη�������ĳһ����Լ���������������ĩβԪ��û�в����Լ
	��������Сֵ�����Բ��ϴ����
	*/
	myint block_threads = GPU::get_max_threads_per_block();// ÿ��block���߳���
	myint threads_needed = n / 2;// n��Ԫ�أ���1�ι�Լ��Ҫn/2���߳�
	myint blocks = threads_needed / block_threads + (threads_needed % block_threads > 0 ? 1 : 0);
	return blocks;
}

// ��ȡ���ڵ���original_num����С��2���ݴη�

myint GPU::ReduceHelper::get_next_power_2_number(myint n) {
	/*
	������2,147,483,647ʱ��next_pow2���ջ�����1,073,741,824��Ȼ����-2,147,483,648���ٳ���2���0��������ѭ��
	�����Ҫ�ж�myint�����ֵ
	size of int: 4��������=2^(4*8)-1=2,147,483,647
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
