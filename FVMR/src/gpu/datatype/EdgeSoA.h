#ifndef _EDGE_SOA_H
#define _EDGE_SOA_H

#include "DefineType.h"
#include "../Env.h"

namespace GPU {
	// structure of array
	struct EdgeSoA {
	public:
		myint num_edge = 0;// 仅可CPU读取
		myint* ID;
		myint* nodes[2];
		myint* setID;// 存储边界类型ID。虽然边界不会超过21亿种，但为了内存对齐，全部用myint
		myint* periodicPair;// 周期边界对应pair的ID。普通边界默认是-1
		myint* elementL;// 左单元ID。任何edge一定有左单元，因此一定不是-1
		myint* elementR;// 右单元ID。边界edge的右单元一般是-1，周期边界除外
		myfloat* length;
		myfloat* distanceOfElements;// 两侧单元中心距离
		myfloat* xy[2];// 边坐标
		myfloat* normal[2];// 边法向量(归一化)

	public:
		// 申请内存 初始化放在后面，因为后面被覆盖，初始化没有意义
		void alloc(myint num) {
			num_edge = num;

			ID = new myint[num];
			nodes[0] = new myint[num];
			nodes[1] = new myint[num];
			setID = new myint[num];
			periodicPair = new myint[num];
			elementL = new myint[num];
			elementR = new myint[num];
			length = new myfloat[num];
			distanceOfElements = new myfloat[num];
			xy[0] = new myfloat[num];
			xy[1] = new myfloat[num];
			normal[0] = new myfloat[num];
			normal[1] = new myfloat[num];
		}

		void free() {
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
		void cuda_alloc(myint num) {
			num_edge = num;

			cudaMalloc((void**)&ID, num * sizeof(myint));
			cudaMalloc((void**)&nodes[0], num * sizeof(myint));
			cudaMalloc((void**)&nodes[1], num * sizeof(myint));
			cudaMalloc((void**)&setID, num * sizeof(myint));
			cudaMalloc((void**)&periodicPair, num * sizeof(myint));
			cudaMalloc((void**)&elementL, num * sizeof(myint));
			cudaMalloc((void**)&elementR, num * sizeof(myint));
			cudaMalloc((void**)&length, num * sizeof(myfloat));
			cudaMalloc((void**)&distanceOfElements, num * sizeof(myfloat));
			cudaMalloc((void**)&xy[0], num * sizeof(myfloat));
			cudaMalloc((void**)&xy[1], num * sizeof(myfloat));
			cudaMalloc((void**)&normal[0], num * sizeof(myfloat));
			cudaMalloc((void**)&normal[1], num * sizeof(myfloat));
		}

		void cuda_free() {
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

		bool has(myint iEdge) {
			return (iEdge >= 0 && iEdge < num_edge);
		}
	};

}

#endif // !_EDGE_SOA_H
