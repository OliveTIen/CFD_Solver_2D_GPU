#ifndef _NODE_SOA_H_
#define _NODE_SOA_H_

#include "DefineType.h"
#include "../Env.h"
// 该类仅提供基本数据类型，初始化、释放、拷贝等操作，不提供任何计算函数，
// 以实现数据和方法的分离。计算功能见
namespace GPU {
	struct NodeSoA {
	public:
		myint num_node = 0;

		myint* ID;
		myint* neighbor_element_count;
		myfloat* xy[2];
	public:
		void alloc(myint _num_node) {
			num_node = _num_node;

			ID = new myint[num_node];
			neighbor_element_count = new myint[num_node];
			xy[0] = new myfloat[num_node];
			xy[1] = new myfloat[num_node];
		}

		void free() {

			delete[] ID;
			delete[] neighbor_element_count;
			delete[] xy[0];
			delete[] xy[1];
		}

		void cuda_alloc(myint _num_node) {
			num_node = _num_node;

			cudaMalloc(&ID, num_node * sizeof(myint));
			cudaMalloc(&neighbor_element_count, num_node * sizeof(myint));
			cudaMalloc(&xy[0], num_node * sizeof(myfloat));
			cudaMalloc(&xy[1], num_node * sizeof(myfloat));
		}

		void cuda_free() {

			cudaFree(ID);
			cudaFree(neighbor_element_count);
			cudaFree(xy[0]);
			cudaFree(xy[1]);
		}

		static void cuda_memcpy(NodeSoA* dst, const NodeSoA* src, cudaMemcpyKind kind);
	};
}

#endif