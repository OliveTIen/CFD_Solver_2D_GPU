#ifndef _NODE_SOA_H_
#define _NODE_SOA_H_

#include "Define.h"
// 该类仅提供基本数据类型，初始化、释放、拷贝等操作，不提供任何计算函数，
// 以实现数据和方法的分离。计算功能见
namespace GPU {
	class NodeSoA {
	public:
		int num_node = 0;

		int* ID;
		REAL* xy[2];
	public:
		void alloc(int _num_node) {
			num_node = _num_node;
			ID = new int[num_node];
			xy[0] = new REAL[num_node];
			xy[1] = new REAL[num_node];
		}

		void free() {
			delete[] ID;
			delete[] xy[0];
			delete[] xy[1];
		}

		void cuda_alloc(int _num_node) {
			num_node = _num_node;
			cudaMalloc(&ID, num_node * sizeof(int));
			cudaMalloc(&xy[0], num_node * sizeof(REAL));
			cudaMalloc(&xy[1], num_node * sizeof(REAL));
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