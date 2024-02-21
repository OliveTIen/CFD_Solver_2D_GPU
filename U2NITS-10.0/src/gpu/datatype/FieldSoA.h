#ifndef _FIELD_SOA_H_
#define _FIELD_SOA_H_

#include "Define.h"

namespace GPU {
	/// <summary>
	/// 场变量 适合node和element
	/// 
	/// </summary>
	class FieldSoA {
	public:
		int num;

		REAL* U[4];
		REAL* Ux[4];
		REAL* Uy[4];
		REAL* Flux[4];

	public:
		void alloc(int _num) {
			num = _num;
			for (int i = 0; i < 4; i++) {
				U[i] = new REAL[num];
				Ux[i] = new REAL[num];
				Uy[i] = new REAL[num];
				Flux[i] = new REAL[num];
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

		void cuda_alloc(int _num) {
			num = _num;
			for (int i = 0; i < 4; i++) {
				cudaMalloc(&U[i], num * sizeof(REAL));
				cudaMalloc(&Ux[i], num * sizeof(REAL));
				cudaMalloc(&Uy[i], num * sizeof(REAL));
				cudaMalloc(&Flux[i], num * sizeof(REAL));
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

		static void cuda_memcpy(FieldSoA* dist, const FieldSoA* src, cudaMemcpyKind kind);
	};

	// 面守恒量 待完成
	class EdgeFieldSoA {
	public:
		int num_edge;

	};

	// 用于输出Tecplot文件的场变量
	class OutputNodeFieldSoA {
	public:
		int num_node;

		REAL* ruvp[4];

	public:
		void alloc(int _num_node) {
			num_node = _num_node;
			for (int i = 0; i < 4; i++) {
				ruvp[i] = new REAL[num_node];
			}
		}
		void free() {
			for (int i = 0; i < 4; i++) {
				delete[] ruvp[i];
			}
		}
		void cuda_alloc(int _num_node) {
			num_node = _num_node;
			for (int i = 0; i < 4; i++) {
				cudaMalloc(&ruvp[i], num_node * sizeof(REAL));
			}
		}
		void cuda_free() {
			for (int i = 0; i < 4; i++) {
				cudaFree(ruvp[i]);
			}
		}
		static void cuda_memcpy(OutputNodeFieldSoA* dist, const OutputNodeFieldSoA* src, cudaMemcpyKind kind);
	};

}


#endif