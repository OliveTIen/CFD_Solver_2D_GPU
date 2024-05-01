#ifndef _FIELD_SOA_H_
#define _FIELD_SOA_H_

#include "DefineType.h"
#include "../Env.h"
#include <vector>

namespace GPU {
	/// <summary>
	/// 场变量 适合node和element
	/// 
	/// </summary>
	struct ElementFieldSoA {
	public:
		myint num;

		myfloat* U[4];
		myfloat* Ux[4];
		myfloat* Uy[4];
		myfloat* Flux[4];

	public:
		void alloc(myint _num) {
			num = _num;

			for (int i = 0; i < 4; i++) {
				U[i] = new myfloat[num];
				Ux[i] = new myfloat[num];
				Uy[i] = new myfloat[num];
				Flux[i] = new myfloat[num];
			}
		}

		void copyfrom(const ElementFieldSoA& src) {
			for (int i = 0; i < 4; i++) {
				for (myint j = 0; j < num; j++) {
					U[i][j] = src.U[i][j];
					Ux[i][j] = src.Ux[i][j];
					Uy[i][j] = src.Uy[i][j];
					Flux[i][j] = src.Flux[i][j];
				}
			}
		}

		void free() {

			for (int i = 0; i < 4; i++) {
				delete[] U[i];
				delete[] Ux[i];
				delete[] Uy[i];
				delete[] Flux[i];
			}
		}

		void cuda_alloc(myint _num) {
			num = _num;


			for (int i = 0; i < 4; i++) {
				cudaMalloc(&U[i], num * sizeof(myfloat));
				cudaMalloc(&Ux[i], num * sizeof(myfloat));
				cudaMalloc(&Uy[i], num * sizeof(myfloat));
				cudaMalloc(&Flux[i], num * sizeof(myfloat));
			}
		}

		void cuda_free() {

			for (int i = 0; i < 4; i++) {
				cudaFree(U[i]);
				cudaFree(Ux[i]);
				cudaFree(Uy[i]);
				cudaFree(Flux[i]);
			}
		}

		static void cuda_memcpy(ElementFieldSoA* dist, const ElementFieldSoA* src, cudaMemcpyKind kind);

		bool has(myint iElement) {
			return (iElement >= 0 && iElement < num);
		}
	};

	// 面守恒量 待完成
	struct EdgeFieldSoA {

		myint num_edge;
		myfloat* Flux[4];

		void alloc(myint _num_edge) {
			num_edge = _num_edge;
			for (int i = 0; i < 4; i++) {
				Flux[i] = new myfloat[num_edge];
			}
		}
		void free() {
			for (int i = 0; i < 4; i++) {
				delete[] Flux[i];
			}
		}
		void cuda_alloc(myint _num_edge) {
			num_edge = _num_edge;
			for (int i = 0; i < 4; i++) {
				cudaMalloc(&Flux[i], num_edge * sizeof(myfloat));
			}
		}
		void cuda_free() {
			for (int i = 0; i < 4; i++) {
				cudaFree(Flux[i]);
			}
		}
		static void cuda_memcpy(EdgeFieldSoA* dist, const EdgeFieldSoA* src, cudaMemcpyKind kind);
	};

	/*
	用于输出Tecplot文件的场变量。
	ruvp采用传统的alloc和free，因为一定会输出到流场
	但其他变量要根据config确定是否输出，因此用动态数组

	初始化时，FieldWriter根据其OutputScheme确定vector大小
	参见FieldWriter::allocDataUsingOutputScheme
	*/
	class OutputNodeFieldSoA {
	public:
		bool b_ruvp_allocated = false;// ruvp是否已经申请内存
		bool b_all_allocated = false;// 待输出的变量数组是否已申请内存
		myint num_node = 0;

		myfloat* ruvp[4]{nullptr, nullptr, nullptr, nullptr};
		std::vector<myfloat> Cp;
		std::vector<myfloat> Ma;


	public:
		// 申请ruvp资源，使用前应判断是否重复申请。由于LogWriter不适合放在头文件，应把判断放在caller中
		void alloc_ruvp(myint _num_node) {
			num_node = _num_node;
			for (int i = 0; i < 4; i++) {
				ruvp[i] = new myfloat[num_node];
			}
			b_ruvp_allocated = true;
		}
		void free_ruvp() {
			for (int i = 0; i < 4; i++) {
				delete[] ruvp[i];
			}
			b_ruvp_allocated = false;
		}
	};

}


#endif