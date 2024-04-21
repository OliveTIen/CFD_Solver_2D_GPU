
#ifndef _BOUNDARY_SET_MAP_H_
#define _BOUNDARY_SET_MAP_H_
#include "DefineType.h"
#include "../Env.h"
namespace GPU {
	struct BoundarySetMap {
	public:
		integer size = 0;

		integer* type = nullptr;

		// 以下两个待完成。存储各个边界有哪些edge。便于遍历边界边，进而积分气动力
		//integer* numOfEdge = nullptr;
		//integer** edge = nullptr;

		// 主要是如果要用指针，传参比较麻烦，既要传num of boundary set，又要传每个set的edge的数量。但是
		// 我不想立刻申请edge数量的内存，想等到复制时再传。
		// 目前的妥协方式是先用BoundaryManager_2D
	public:
		void alloc(integer _size) {
			size = _size;

			type = new integer[size];
			//numOfEdge = new integer[size];
			//for (integer i = 0; i < size; i++) {
			//	numOfEdge[i]=new integer
			//}
		}
		void free() {

			delete[] type;
			//delete[] numOfEdge;
		}
		void cuda_alloc(integer _size) {
			size = _size;

			cudaMalloc(&type, size * sizeof(integer));
			//cudaMalloc(&numOfEdge, size * sizeof(integer));
		}

		void cuda_free() {

			cudaFree(type);
			//cudaFree(numOfEdge);
		}

		static void cuda_memcpy(BoundarySetMap* dist, const BoundarySetMap* src, cudaMemcpyKind kind);

	};
}

#endif // !_BOUNDARY_SET_MAP_H_
