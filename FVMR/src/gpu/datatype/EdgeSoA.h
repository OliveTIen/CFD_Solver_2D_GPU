#ifndef _EDGE_SOA_H
#define _EDGE_SOA_H

#include "DefineType.h"
#include "../Env.h"

namespace GPU {
	// structure of array
	struct EdgeSoA {
	public:
		myint num_edge = 0;// ����CPU��ȡ
		myint* ID;
		myint* nodes[2];
		myint* setID;// �洢�߽�����ID����Ȼ�߽粻�ᳬ��21���֣���Ϊ���ڴ���룬ȫ����myint
		myint* periodicPair;// ���ڱ߽��Ӧpair��ID����ͨ�߽�Ĭ����-1
		myint* elementL;// ��ԪID���κ�edgeһ������Ԫ�����һ������-1
		myint* elementR;// �ҵ�ԪID���߽�edge���ҵ�Ԫһ����-1�����ڱ߽����
		myfloat* length;
		myfloat* distanceOfElements;// ���൥Ԫ���ľ���
		myfloat* xy[2];// ������
		myfloat* normal[2];// �߷�����(��һ��)

	public:
		// �����ڴ� ��ʼ�����ں��棬��Ϊ���汻���ǣ���ʼ��û������
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
		// cuda �����ڴ�
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
