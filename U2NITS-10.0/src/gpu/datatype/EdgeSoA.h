#ifndef _EDGE_SOA_H
#define _EDGE_SOA_H

#include "Define.h"

namespace GPU {
	// structure of array
	struct EdgeSoA {
	public:
		int num_edge = 0;// ����CPU��ȡ
		int* _num_edge_;// ͨ�������ڴ�õ��ģ�����GPU��ȡ��������cuda_memcpy()����Ϊ��ʼ��ʱ��ȷ��ֵ
		int* ID;
		int* nodes[2];
		int* setID;
		int* periodicPair;// ���ڱ߽��Ӧpair��ID����ͨ�߽�Ĭ����-1
		int* elementL;// ��ԪID���κ�edgeһ������Ԫ�����һ������-1
		int* elementR;// �ҵ�ԪID���߽�edge���ҵ�Ԫһ����-1�����ڱ߽����
		REAL* length;
		REAL* distanceOfElements;// ���൥Ԫ���ľ���
		REAL* xy[2];// ������
		REAL* normal[2];// �߷�����(��һ��)

	public:
		// �����ڴ� ��ʼ�����ں��棬��Ϊ���汻���ǣ���ʼ��û������
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
		// cuda �����ڴ�
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
