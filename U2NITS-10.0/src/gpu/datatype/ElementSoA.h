#ifndef _ELEMENT_SOA_H_
#define _ELEMENT_SOA_H_

#include "Define.h"
/*
ElementSoA�����������ݣ�һ���ǻ�����������ݣ���ID���ڵ㡢�ߡ��ھӵȣ�
��һ���Ǳ仯�����ݣ���U��Ux��Uy��Flux�ȡ�
�仯��������Ҫ��host��device֮��ͬ���������Ҫ�õ�cudaMemcpy��


*/

namespace GPU {
	class ElementSoA {
	public:
		int num_element;

		int* ID;
		int* nodes[4];
		int* edges[4];// faces
		int* neighbors[4];// �ھӵ�Ԫ����edges˳�����У�-1��ʾ���ھӣ����ڱ߽���Ȼ��-1
		REAL* xy[2];
		REAL* volume;
		//REAL* y;

		//REAL* U[4];
		//REAL* Ux[4];
		//REAL* Uy[4];
		//REAL* Flux[4];
		
	public:
		void alloc(int _num_element) {
			num_element = _num_element;
			ID = new int[num_element];
			for (int i = 0; i < 4; i++) {
				nodes[i] = new int[num_element];
				edges[i] = new int[num_element];
				neighbors[i] = new int[num_element];
				//U[i] = new REAL[num_element];
				//Ux[i] = new REAL[num_element];
				//Uy[i] = new REAL[num_element];
				//Flux[i] = new REAL[num_element];
			}
			xy[0] = new REAL[num_element];
			xy[1] = new REAL[num_element];
			volume = new REAL[num_element];
		}

		void free() {
			delete[] ID;
			for (int i = 0; i < 4; i++) {
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

		void cuda_alloc(int _num_element) {
			num_element = _num_element;
			cudaMalloc(&ID, num_element * sizeof(int));
			for (int i = 0; i < 4; i++) {
				cudaMalloc(&nodes[i], num_element * sizeof(int));
				cudaMalloc(&edges[i], num_element * sizeof(int));
				cudaMalloc(&neighbors[i], num_element * sizeof(int));
				//cudaMalloc(&U[i], num_element * sizeof(REAL));
				//cudaMalloc(&Ux[i], num_element * sizeof(REAL));
				//cudaMalloc(&Uy[i], num_element * sizeof(REAL));
				//cudaMalloc(&Flux[i], num_element * sizeof(REAL));
			}
			cudaMalloc(&xy[0], num_element * sizeof(REAL));
			cudaMalloc(&xy[1], num_element * sizeof(REAL));
			cudaMalloc(&volume, num_element * sizeof(REAL));
		}

		void cuda_free() {
			cudaFree(ID);
			for (int i = 0; i < 4; i++) {
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

	};
}

/*

public:
	//�ṹ��Ϣ
	int ID = -1;
	int node_num = 3;//��Ԫ����
	int nodes[4]{-1,-1,-1,-1};//node ID ֮������ID������ָ�룬����Ϊ��ȡ�ļ���ֱ�Ӷ���ElementID-NodeID-NodeID-NodeID
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
	int GPUID = -1;

*/

#endif