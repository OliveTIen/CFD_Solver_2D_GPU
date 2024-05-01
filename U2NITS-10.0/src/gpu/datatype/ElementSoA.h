#ifndef _ELEMENT_SOA_H_
#define _ELEMENT_SOA_H_

#include "DefineType.h"
#include "../Env.h"
/*
ElementSoA�����������ݣ�һ���ǻ�����������ݣ���ID���ڵ㡢�ߡ��ھӵȣ�
��һ���Ǳ仯�����ݣ���U��Ux��Uy��Flux�ȡ�
�仯��������Ҫ��host��device֮��ͬ���������Ҫ�õ�cudaMemcpy��


*/

namespace GPU {
	struct ElementSoA {
	public:
		myint num_element = 0;

		myint* ID;
		myint* nodes[4];
		myint* edges[4];// faces�������ε�Ԫ����4����ȡ-1
		myint* neighbors[4];// �ھӵ�Ԫ����edges˳�����У�-1��ʾ���ھӣ����ڱ߽���Ȼ��-1
		myfloat* xy[2];
		myfloat* volume;
		//myfloat* y;

		//myfloat* U[4];
		//myfloat* Ux[4];
		//myfloat* Uy[4];
		//myfloat* Flux[4];

	public:
		void alloc(myint _num_element) {
			num_element = _num_element;
			ID = new myint[num_element];
			for (myint i = 0; i < 4; i++) {
				nodes[i] = new myint[num_element];
				edges[i] = new myint[num_element];
				neighbors[i] = new myint[num_element];
				//U[i] = new myfloat[num_element];
				//Ux[i] = new myfloat[num_element];
				//Uy[i] = new myfloat[num_element];
				//Flux[i] = new myfloat[num_element];
			}
			xy[0] = new myfloat[num_element];
			xy[1] = new myfloat[num_element];
			volume = new myfloat[num_element];
		}

		void free() {
			delete[] ID;
			for (myint i = 0; i < 4; i++) {
				delete[] nodes[i];
				delete[] edges[i];
				delete[] neighbors[i];
				//delete[] U[i];
				//delete[] Ux[i];
				//delete[] Uy[i];
				//delete[] Flux[i];
			}
			delete[] xy[0];
			delete[] xy[1];
			delete[] volume;
		}

		void cuda_alloc(myint _num_element) {
			num_element = _num_element;

			cudaMalloc(&ID, num_element * sizeof(myint));
			for (myint i = 0; i < 4; i++) {
				cudaMalloc(&nodes[i], num_element * sizeof(myint));
				cudaMalloc(&edges[i], num_element * sizeof(myint));
				cudaMalloc(&neighbors[i], num_element * sizeof(myint));
				//cudaMalloc(&U[i], num_element * sizeof(myfloat));
				//cudaMalloc(&Ux[i], num_element * sizeof(myfloat));
				//cudaMalloc(&Uy[i], num_element * sizeof(myfloat));
				//cudaMalloc(&Flux[i], num_element * sizeof(myfloat));
			}
			cudaMalloc(&xy[0], num_element * sizeof(myfloat));
			cudaMalloc(&xy[1], num_element * sizeof(myfloat));
			cudaMalloc(&volume, num_element * sizeof(myfloat));
		}

		void cuda_free() {
			cudaFree(ID);
			for (myint i = 0; i < 4; i++) {
				cudaFree(nodes[i]);
				cudaFree(edges[i]);
				cudaFree(neighbors[i]);
				//cudaFree(U[i]);
				//cudaFree(Ux[i]);
				//cudaFree(Uy[i]);
				//cudaFree(Flux[i]);
			}
			cudaFree(xy[0]);
			cudaFree(xy[1]);
			cudaFree(volume);
		}

		static void cuda_memcpy(ElementSoA* dist, const ElementSoA* src, cudaMemcpyKind kind);

		bool has(myint iElement) {
			return (iElement >= 0 && iElement < num_element);
		}
	};
}

/*

public:
	//�ṹ��Ϣ
	myint ID = -1;
	myint node_num = 3;//��Ԫ����
	myint nodes[4]{-1,-1,-1,-1};//node ID ֮������ID������ָ�룬����Ϊ��ȡ�ļ���ֱ�Ӷ���ElementID-NodeID-NodeID-NodeID
	Edge_2D* pEdges[4]{};//�Ѿ��ڶ�ȡ�ļ�ʱ��ʼ����������
	//������Ϣ
	double x = 0;//��Ԫ�������ꡣʹ��ǰ���calxy!
	double y = 0;//��Ԫ��������
	double U[4] = { 1,0,0,1.1e5 / 0.4 };//�غ�����,��u,��v,��E
	double Ux[4] = { 0,0,0,0 };//
	double Uy[4] = { 0,0,0,0 };//
	double Flux[4]{};//4���غ�������ֵͨ����ÿ�μӼ�ǰҪ����
	//������Ϣ
	double deltaeig;//Roe��ʽ����ͨ��ʱ�õ���ÿһ��������

	//static Eigen::MatrixXi si_ti;//����Ĳ�������
	//static Eigen::MatrixXd GaussPointMatrix;//����Ĳ�������
	myint GPUID = -1;

*/

#endif