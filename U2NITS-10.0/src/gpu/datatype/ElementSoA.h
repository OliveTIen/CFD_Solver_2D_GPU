#ifndef _ELEMENT_SOA_H_
#define _ELEMENT_SOA_H_

#include "Define.h"
/*
ElementSoA中有两种数据，一种是基本不变的数据，如ID、节点、边、邻居等，
另一种是变化的数据，如U、Ux、Uy、Flux等。
变化的数据需要在host和device之间同步，因此需要用到cudaMemcpy。


*/

namespace GPU {
	struct ElementSoA {
	public:
		integer num_element = 0;
		integer* _num_element_;// GPU使用

		integer* ID;
		integer* nodes[4];
		integer* edges[4];// faces
		integer* neighbors[4];// 邻居单元，按edges顺序排列；-1表示无邻居；周期边界依然是-1
		REAL* xy[2];
		REAL* volume;
		//REAL* y;

		//REAL* U[4];
		//REAL* Ux[4];
		//REAL* Uy[4];
		//REAL* Flux[4];
		
	public:
		void alloc(integer _num_element) {
			num_element = _num_element;
			_num_element_ = new integer;
			*_num_element_ = num_element;
			ID = new integer[num_element];
			for (integer i = 0; i < 4; i++) {
				nodes[i] = new integer[num_element];
				edges[i] = new integer[num_element];
				neighbors[i] = new integer[num_element];
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
			delete _num_element_;
			delete[] ID;
			for (integer i = 0; i < 4; i++) {
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

		void cuda_alloc(integer _num_element) {
			num_element = _num_element;
			cudaMalloc(&_num_element_, 1 * sizeof(integer));
			cudaMemcpy(_num_element_, &num_element, 1 * sizeof(integer), ::cudaMemcpyHostToDevice);

			cudaMalloc(&ID, num_element * sizeof(integer));
			for (integer i = 0; i < 4; i++) {
				cudaMalloc(&nodes[i], num_element * sizeof(integer));
				cudaMalloc(&edges[i], num_element * sizeof(integer));
				cudaMalloc(&neighbors[i], num_element * sizeof(integer));
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
			cudaFree(_num_element_);
			cudaFree(ID);
			for (integer i = 0; i < 4; i++) {
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

		bool has(integer iElement) {
			return (iElement >= 0 && iElement < num_element);
		}
	};
}

/*

public:
	//结构信息
	integer ID = -1;
	integer node_num = 3;//单元类型
	integer nodes[4]{-1,-1,-1,-1};//node ID 之所以用ID而不是指针，是因为读取文件是直接读的ElementID-NodeID-NodeID-NodeID
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
	integer GPUID = -1;

*/

#endif