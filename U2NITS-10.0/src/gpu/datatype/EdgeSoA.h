#ifndef _EDGE_SOA_H
#define _EDGE_SOA_H

#include "Define.h"

namespace GPU {
	// structure of array
	struct EdgeSoA {
	public:
		int num_edge = 0;// 仅可CPU读取
		int* _num_edge_;// 通过申请内存得到的，可以GPU读取。不参与cuda_memcpy()，因为初始化时已确定值
		int* ID;
		int* nodes[2];
		int* setID;
		int* periodicPair;// 周期边界对应pair的ID。普通边界默认是-1
		int* elementL;// 左单元ID。任何edge一定有左单元，因此一定不是-1
		int* elementR;// 右单元ID。边界edge的右单元一般是-1，周期边界除外
		REAL* length;
		REAL* distanceOfElements;// 两侧单元中心距离
		REAL* xy[2];// 边坐标
		REAL* normal[2];// 边法向量(归一化)

	public:
		// 申请内存 初始化放在后面，因为后面被覆盖，初始化没有意义
		void alloc(int num) {
			num_edge = num;
			_num_edge_ = new int;
			*_num_edge_ = num_edge;

			ID = new int[num];
			nodes[0] = new int[num];
			nodes[1] = new int[num];
			setID = new int[num];
			periodicPair = new int[num];
			elementL = new int[num];
			elementR = new int[num];
			length = new REAL[num];
			distanceOfElements = new REAL[num];	
			xy[0] = new REAL[num];
			xy[1] = new REAL[num];
			normal[0] = new REAL[num];
			normal[1] = new REAL[num];
		}

		void free() {
			delete _num_edge_;
			delete[] ID;
			delete[] nodes[0];
			delete[] nodes[1];
			delete[] setID;
			delete[] periodicPair;
			delete[] elementL;
			delete[] elementR;
			delete[] length;
			delete[] distanceOfElements;
			delete[] xy[0];
			delete[] xy[1];
			delete[] normal[0];
			delete[] normal[1];
			num_edge = 0;
		}
		// cuda 申请内存
		void cuda_alloc(int num) {
			num_edge = num;
			cudaMalloc((void**)&_num_edge_, sizeof(int));
			cudaMemcpy(_num_edge_, &num_edge, sizeof(int), ::cudaMemcpyHostToDevice);

			cudaMalloc((void**)&ID, num * sizeof(int));
			cudaMalloc((void**)&nodes[0], num * sizeof(int));
			cudaMalloc((void**)&nodes[1], num * sizeof(int));
			cudaMalloc((void**)&setID, num * sizeof(int));
			cudaMalloc((void**)&periodicPair, num * sizeof(int));
			cudaMalloc((void**)&elementL, num * sizeof(int));
			cudaMalloc((void**)&elementR, num * sizeof(int));
			cudaMalloc((void**)&length, num * sizeof(REAL));
			cudaMalloc((void**)&distanceOfElements, num * sizeof(REAL));
			cudaMalloc((void**)&xy[0], num * sizeof(REAL));
			cudaMalloc((void**)&xy[1], num * sizeof(REAL));
			cudaMalloc((void**)&normal[0], num * sizeof(REAL));
			cudaMalloc((void**)&normal[1], num * sizeof(REAL));
		}

		void cuda_free() {
			cudaFree(_num_edge_);
			cudaFree(ID);
			cudaFree(nodes[0]);
			cudaFree(nodes[1]);
			cudaFree(setID);
			cudaFree(periodicPair);
			cudaFree(elementL);
			cudaFree(elementR);
			cudaFree(length);
			cudaFree(distanceOfElements);
			cudaFree(xy[0]);
			cudaFree(xy[1]);
			cudaFree(normal[0]);
			cudaFree(normal[1]);
			num_edge = 0;
		}

		static void cuda_memcpy(EdgeSoA* dist, const EdgeSoA* src, cudaMemcpyKind kind);

		bool has(integer iEdge) {
			return (iEdge >= 0 && iEdge < num_edge);
		}
	};

}

#endif // !_EDGE_SOA_H
