#ifndef _NODE_SOA_H_
#define _NODE_SOA_H_

#include "Define.h"
// ������ṩ�����������ͣ���ʼ�����ͷš������Ȳ��������ṩ�κμ��㺯����
// ��ʵ�����ݺͷ����ķ��롣���㹦�ܼ�
namespace GPU {
	class NodeSoA {
	public:
		int num_node = 0;

		int* ID;
		REAL* x;
		REAL* y;
	public:
		void alloc(int _num_node) {
			num_node = _num_node;
			ID = new int[num_node];
			x = new REAL[num_node];
			y = new REAL[num_node];
		}

		void free() {
			delete[] ID;
			delete[] x;
			delete[] y;
		}

		void cuda_alloc(int _num_node) {
			num_node = _num_node;
			cudaMalloc(&ID, num_node * sizeof(int));
			cudaMalloc(&x, num_node * sizeof(REAL));
			cudaMalloc(&y, num_node * sizeof(REAL));
		}

		void cuda_free() {
			cudaFree(ID);
			cudaFree(x);
			cudaFree(y);
		}

		static void cuda_memcpy(NodeSoA* dist, const NodeSoA* src, cudaMemcpyKind kind);
	};
}

#endif