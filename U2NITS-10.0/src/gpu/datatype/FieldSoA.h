#ifndef _FIELD_SOA_H_
#define _FIELD_SOA_H_

#include "DefineType.h"
#include "../Env.h"

namespace GPU {
	/// <summary>
	/// 场变量 适合node和element
	/// 
	/// </summary>
	struct FieldSoA {
	public:
		int num;

		real* U[4];
		real* Ux[4];
		real* Uy[4];
		real* Flux[4];

	public:
		void alloc(int _num) {
			num = _num;

			for (int i = 0; i < 4; i++) {
				U[i] = new real[num];
				Ux[i] = new real[num];
				Uy[i] = new real[num];
				Flux[i] = new real[num];
			}
		}

		void copyfrom(const FieldSoA& src) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < num; j++) {
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

		void cuda_alloc(int _num) {
			num = _num;


			for (int i = 0; i < 4; i++) {
				cudaMalloc(&U[i], num * sizeof(real));
				cudaMalloc(&Ux[i], num * sizeof(real));
				cudaMalloc(&Uy[i], num * sizeof(real));
				cudaMalloc(&Flux[i], num * sizeof(real));
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
		real* numeralFlux[4];

	public:
		void alloc(int _num_edge) {
			num_edge = _num_edge;
			for (int i = 0; i < 4; i++) {
				numeralFlux[i] = new real[num_edge];
			}
		}
		void free() {
			for (int i = 0; i < 4; i++) {
				delete[] numeralFlux[i];
			}
		}
		void cuda_alloc(int _num_edge) {
			num_edge = _num_edge;
			for (int i = 0; i < 4; i++) {
				cudaMalloc(&numeralFlux[i], num_edge * sizeof(real));
			}
		}
		void cuda_free() {
			for (int i = 0; i < 4; i++) {
				cudaFree(numeralFlux[i]);
			}
		}
		static void cuda_memcpy(EdgeFieldSoA* dist, const EdgeFieldSoA* src, cudaMemcpyKind kind);
	};

	// 用于输出Tecplot文件的场变量
	class OutputNodeFieldSoA {
	public:
		int num_node;

		real* ruvp[4];

	public:
		void alloc(int _num_node) {
			num_node = _num_node;
			for (int i = 0; i < 4; i++) {
				ruvp[i] = new real[num_node];
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
				cudaMalloc(&ruvp[i], num_node * sizeof(real));
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