#ifndef _NODE_SOA_H_
#define _NODE_SOA_H_

#include "DefineType.h"
#include "../Env.h"
// ������ṩ�����������ͣ���ʼ�����ͷš������Ȳ��������ṩ�κμ��㺯����
// ��ʵ�����ݺͷ����ķ��롣���㹦�ܼ�
namespace GPU {
	struct NodeSoA {
	public:
		int num_node = 0;

		int* ID;
		myfloat* xy[2];
	public:
		void alloc(int _num_node) {
			num_node = _num_node;

			ID = new int[num_node];
			xy[0] = new myfloat[num_node];
			xy[1] = new myfloat[num_node];
		}

		void free() {

			delete[] ID;
			delete[] xy[0];
			delete[] xy[1];
		}

		void cuda_alloc(int _num_node) {
			num_node = _num_node;

			cudaMalloc(&ID, num_node * sizeof(int));
			cudaMalloc(&xy[0], num_node * sizeof(myfloat));
			cudaMalloc(&xy[1], num_node * sizeof(myfloat));
		}

		void cuda_free() {

			cudaFree(ID);
			cudaFree(xy[0]);
			cudaFree(xy[1]);
		}

		static void cuda_memcpy(NodeSoA* dist, const NodeSoA* src, cudaMemcpyKind kind);
	};
}

#endif