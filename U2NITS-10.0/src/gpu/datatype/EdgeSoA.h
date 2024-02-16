#ifndef _EDGE_SOA_H
#define _EDGE_SOA_H

#include "Define.h"

namespace GPU {
	// structure of array
	class EdgeSoA {
	public:
		int num_edge = 0;
		int* ID;
		int* nodes[2];
		int* setID;
		int* elementL;
		int* elementR;
		REAL* length;
		REAL* distanceOfElements;// 两侧单元中心距离

	public:
		// 申请内存 初始化放在后面，因为后面被覆盖，初始化没有意义
		void alloc(int num) {
			num_edge = num;
			ID = new int[num];
			nodes[0] = new int[num];
			nodes[1] = new int[num];
			setID = new int[num];
			elementL = new int[num];
			elementR = new int[num];
			length = new REAL[num];
			distanceOfElements = new REAL[num];	
		}

		void free() {
			delete[] ID;
			delete[] nodes[0];
			delete[] nodes[1];
			delete[] setID;
			delete[] elementL;
			delete[] elementR;
			delete[] length;
			delete[] distanceOfElements;
			num_edge = 0;
		}
		// cuda 申请内存
		void cuda_alloc(int num) {
			num_edge = num;
			cudaMalloc((void**)&ID, num * sizeof(int));
			cudaMalloc((void**)&nodes[0], num * sizeof(int));
			cudaMalloc((void**)&nodes[1], num * sizeof(int));
			cudaMalloc((void**)&setID, num * sizeof(int));
			cudaMalloc((void**)&elementL, num * sizeof(int));
			cudaMalloc((void**)&elementR, num * sizeof(int));
			cudaMalloc((void**)&length, num * sizeof(REAL));
			cudaMalloc((void**)&distanceOfElements, num * sizeof(REAL));
		}

		void cuda_free() {
			cudaFree(ID);
			cudaFree(nodes[0]);
			cudaFree(nodes[1]);
			cudaFree(setID);
			cudaFree(elementL);
			cudaFree(elementR);
			cudaFree(length);
			cudaFree(distanceOfElements);
			num_edge = 0;
		}

		static void cuda_memcpy(EdgeSoA* dist, const EdgeSoA* src, cudaMemcpyKind kind);
	};

}

#endif // !_EDGE_H
