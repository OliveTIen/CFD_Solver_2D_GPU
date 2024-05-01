#ifndef _ELEMENT_SOA_H_
#define _ELEMENT_SOA_H_

#include "DefineType.h"
#include "../Env.h"
/*
ElementSoA中有两种数据，一种是基本不变的数据，如ID、节点、边、邻居等，
另一种是变化的数据，如U、Ux、Uy、Flux等。
变化的数据需要在host和device之间同步，因此需要用到cudaMemcpy。


*/

namespace GPU {
	struct ElementSoA {
	public:
		myint num_element = 0;

		myint* ID;
		myint* nodes[4];
		myint* edges[4];// faces。三角形单元，第4条边取-1
		myint* neighbors[4];// 邻居单元，按edges顺序排列；-1表示无邻居；周期边界依然是-1
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
	//结构信息
	myint ID = -1;
	myint node_num = 3;//单元类型
	myint nodes[4]{-1,-1,-1,-1};//node ID 之所以用ID而不是指针，是因为读取文件是直接读的ElementID-NodeID-NodeID-NodeID
	Edge_2D* pEdges[4]{};//已经在读取文件时初始化，放心用
	//数据信息
	double x = 0;//单元中心坐标。使用前务必calxy!
	double y = 0;//单元中心坐标
	double U[4] = { 1,0,0,1.1e5 / 0.4 };//守恒量ρ,ρu,ρv,ρE
	double Ux[4] = { 0,0,0,0 };//
	double Uy[4] = { 0,0,0,0 };//
	double Flux[4]{};//4个守恒量的数值通量。每次加减前要清零
	//附加信息
	double deltaeig;//Roe格式计算通量时用到。每一轮需清零

	//static Eigen::MatrixXi si_ti;//顶点的参数坐标
	//static Eigen::MatrixXd GaussPointMatrix;//顶点的参数坐标
	myint GPUID = -1;

*/

#endif