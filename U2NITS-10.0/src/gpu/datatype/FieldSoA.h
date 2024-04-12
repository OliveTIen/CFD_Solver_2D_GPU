#ifndef _FIELD_SOA_H_
#define _FIELD_SOA_H_

#include "Define.h"

namespace GPU {
	/// <summary>
	/// ������ �ʺ�node��element
	/// 
	/// </summary>
	struct FieldSoA {
	public:
		int num;
		int* _num_;// GPUʹ��

		REAL* U[4];
		REAL* Ux[4];
		REAL* Uy[4];
		REAL* Flux[4];

	public:
		void alloc(int _num) {
			num = _num;
			_num_ = new int;
			*_num_ = num;

			for (int i = 0; i < 4; i++) {
				U[i] = new REAL[num];
				Ux[i] = new REAL[num];
				Uy[i] = new REAL[num];
				Flux[i] = new REAL[num];
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
			delete _num_;

			for (int i = 0; i < 4; i++) {
				delete[] U[i];
				delete[] Ux[i];
				delete[] Uy[i];
				delete[] Flux[i];
			}
		}

		void cuda_alloc(int _num) {
			num = _num;
			cudaMalloc(&_num_, 1 * sizeof(int));
			cudaMemcpy(_num_, &num, 1 * sizeof(int), ::cudaMemcpyHostToDevice);


			for (int i = 0; i < 4; i++) {
				cudaMalloc(&U[i], num * sizeof(REAL));
				cudaMalloc(&Ux[i], num * sizeof(REAL));
				cudaMalloc(&Uy[i], num * sizeof(REAL));
				cudaMalloc(&Flux[i], num * sizeof(REAL));
			}
		}

		void cuda_free() {
			cudaFree(_num_);

			for (int i = 0; i < 4; i++) {
				cudaFree(U[i]);
				cudaFree(Ux[i]);
				cudaFree(Uy[i]);
				cudaFree(Flux[i]);
			}
		}

		static void cuda_memcpy(FieldSoA* dist, const FieldSoA* src, cudaMemcpyKind kind);
	};

	// ���غ��� �����
	class EdgeFieldSoA {
	public:
		int num_edge;
		REAL* numeralFlux[4];

	public:
		void alloc(int _num_edge) {
			num_edge = _num_edge;
			for (int i = 0; i < 4; i++) {
				numeralFlux[i] = new REAL[num_edge];
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
				cudaMalloc(&numeralFlux[i], num_edge * sizeof(REAL));
			}
		}
		void cuda_free() {
			for (int i = 0; i < 4; i++) {
				cudaFree(numeralFlux[i]);
			}
		}
		static void cuda_memcpy(EdgeFieldSoA* dist, const EdgeFieldSoA* src, cudaMemcpyKind kind);
	};

	// �������Tecplot�ļ��ĳ�����
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