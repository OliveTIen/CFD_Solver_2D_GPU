
#ifndef _BOUNDARY_SET_MAP_H_
#define _BOUNDARY_SET_MAP_H_
#include "DefineType.h"
#include "../Env.h"
namespace GPU {
	/*
	存储从边界ID到边界类型type的映射。不存储具体有哪些edge
	*/
	struct BoundarySetMap {
	public:
		myint size = 0;

		myint* type = nullptr;

		// 以下两个待完成。存储各个边界有哪些edge。便于遍历边界边，进而积分气动力
		//myint* numOfEdge = nullptr;
		//myint** edge = nullptr;

		// 主要是如果要用指针，传参比较麻烦，既要传num of boundary set，又要传每个set的edge的数量。但是
		// 我不想立刻申请edge数量的内存，想等到复制时再传。
		// 目前的妥协方式是先用BoundaryManager_2D

		/*
		需求：
		BoundarySetMap存储std::vector数组


		问题：如果用std::vector，计算通量时，已知edge，如何获取所属边界类型？
		答案：已经在edge_host中存储了setID，
			先用查询edge_host获取setID：
			int setID = edge_host.setID[iedge];
			再用setID获取边界类型：
			bType = f->boundaryManager.boundaries[setID - 1].type;
		
		问题：如果用std::vector，GPU如何传参？
		答案：GPU中可以用thrust::vector
			先不管GPU，只用CPU
		*/
	public:
		void alloc(myint _size) {
			size = _size;

			type = new myint[size];
			//numOfEdge = new myint[size];
			//for (myint i = 0; i < size; i++) {
			//	numOfEdge[i]=new myint
			//}
		}
		void free() {

			delete[] type;
			//delete[] numOfEdge;
		}
		void cuda_alloc(myint _size) {
			size = _size;

			cudaMalloc(&type, size * sizeof(myint));
			//cudaMalloc(&numOfEdge, size * sizeof(myint));
		}

		void cuda_free() {

			cudaFree(type);
			//cudaFree(numOfEdge);
		}

		static void cuda_memcpy(BoundarySetMap* dist, const BoundarySetMap* src, cudaMemcpyKind kind);

	};
}

#endif // !_BOUNDARY_SET_MAP_H_
