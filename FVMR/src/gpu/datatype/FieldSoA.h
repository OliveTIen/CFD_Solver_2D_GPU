#ifndef _FIELD_SOA_H_
#define _FIELD_SOA_H_

#include "DefineType.h"
#include "../Env.h"
#include <vector>
#include "../GPUGlobalFunction.h"

namespace GPU {
	struct ReduceHelper {
	public:
		// �����Լ������dev_output���鳤��
		static myint get_dev_output_length(myint n);

		// ��ȡ���ڵ���n����С��2���ݴη�
		static myint get_next_power_2_number(myint n);
	};


	struct ElementFieldSoA {
	public:
		myint num;

		myfloat* U[4];
		myfloat* Ux[4];
		myfloat* Uy[4];
		myfloat* Flux[4];

		// 2024-05-16 ���calculateDtGPU��ruvp��host��device������
		myfloat* ruvp[4];
		myfloat* Uold[4];

	public:
		void alloc(myint _num) {
			num = _num;

			for (int i = 0; i < 4; i++) {
				U[i] = new myfloat[num];
				Ux[i] = new myfloat[num];
				Uy[i] = new myfloat[num];
				Flux[i] = new myfloat[num];
				ruvp[i] = new myfloat[num];
				Uold[i] = new myfloat[num];
			}
		}

		void copyfrom(const ElementFieldSoA& src) {
			for (int i = 0; i < 4; i++) {
				for (myint j = 0; j < num; j++) {
					U[i][j] = src.U[i][j];
					Ux[i][j] = src.Ux[i][j];
					Uy[i][j] = src.Uy[i][j];
					Flux[i][j] = src.Flux[i][j];
					ruvp[i][j] = src.Flux[i][j];
					Uold[i][j] = src.Flux[i][j];
				}
			}
		}

		void free() {
			for (int i = 0; i < 4; i++) {
				delete[] U[i];
				delete[] Ux[i];
				delete[] Uy[i];
				delete[] Flux[i];
				delete[] ruvp[i];
				delete[] Uold[i];
			}
		}

		void cuda_alloc(myint _num) {
			num = _num;

			for (int i = 0; i < 4; i++) {
				cudaMalloc(&U[i], num * sizeof(myfloat));
				cudaMalloc(&Ux[i], num * sizeof(myfloat));
				cudaMalloc(&Uy[i], num * sizeof(myfloat));
				cudaMalloc(&Flux[i], num * sizeof(myfloat));
				cudaMalloc(&ruvp[i], num * sizeof(myfloat));
				cudaMalloc(&Uold[i], num * sizeof(myfloat));
			}
		}

		void cuda_free() {

			for (int i = 0; i < 4; i++) {
				cudaFree(U[i]);
				cudaFree(Ux[i]);
				cudaFree(Uy[i]);
				cudaFree(Flux[i]);
				cudaFree(ruvp[i]);
				cudaFree(Uold[i]);
			}
		}

		static void cuda_memcpy(ElementFieldSoA* dst, const ElementFieldSoA* src, cudaMemcpyKind kind);

		bool has(myint iElement) {
			return (iElement >= 0 && iElement < num);
		}
	};

	// ���غ��� �����
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
		static void cuda_memcpy(EdgeFieldSoA* dst, const EdgeFieldSoA* src, cudaMemcpyKind kind);
	};

	/*
	�������Tecplot�ļ��ĳ�������
	ruvp���ô�ͳ��alloc��free����Ϊһ�������������
	����������Ҫ����configȷ���Ƿ����������ö�̬����

	��ʼ��ʱ��FieldWriter������OutputSchemeȷ��vector��С
	�μ�FieldWriter::allocDataUsingOutputScheme
	*/
	struct NodeFieldSoA {
	public:
		myint num_node;
		myfloat* ruvp[4];

		void alloc(myint _num_node) {
			num_node = _num_node;
			for (int i = 0; i < 4; i++) {
				ruvp[i] = new myfloat[num_node];
			}
		}
		void free() {
			for (int i = 0; i < 4; i++) {
				delete[] ruvp[i];
			}
		}
		void cuda_alloc(myint _num_node) {
			num_node = _num_node;
			for (int i = 0; i < 4; i++) {
				cudaMalloc(&ruvp[i], num_node * sizeof(myfloat));
			}
		}
		void cuda_free() {
			for (int i = 0; i < 4; i++) {
				cudaFree(ruvp[i]);
			}
		}
		static void cuda_memcpy(NodeFieldSoA* dst, const NodeFieldSoA* src, cudaMemcpyKind kind);
	};

	// a simple device data. Doesn't have menber cuda_memcpy()
	struct DataSoA {
	public:
		myint length;
		myfloat* data;

		void cuda_alloc(myint _length) {
			length = _length;
			cudaMalloc(&data, length * sizeof(myfloat));
		}
		void cuda_free() {
			cudaFree(data);
		}
	};


	struct OutputNodeFieldSoA {
	public:
		bool b_ruvp_allocated = false;// ruvp�Ƿ��Ѿ������ڴ�
		bool b_all_allocated = false;// ������ı��������Ƿ��������ڴ�
		myint num_node = 0;

		myfloat* ruvp[4]{nullptr, nullptr, nullptr, nullptr};
		std::vector<myfloat> Cp;
		std::vector<myfloat> Ma;


	public:
		// ����ruvp��Դ��ʹ��ǰӦ�ж��Ƿ��ظ����롣����LogWriter���ʺϷ���ͷ�ļ���Ӧ���жϷ���caller��
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

	struct ElementFieldVariable_dt {
		myint num_element = 0;// hostģʽ��alphaC�ĳ���
		myint num_reduce = 0;// ��Լ�����鳤�ȣ���deviceģʽ��alphaC�ĳ��ȣ�ȡ����Ĵ���num_element��2����
		myfloat* alphaC = nullptr;// ����ʱ�������ȴ洢alphaC���ִ洢dt
		myfloat* dev_output = nullptr;// ���ڹ�Լ��������ʱ����������deviceģʽ�������ڴ棬hostģʽ��ֹʹ��

		void alloc(myint _num_element) {
			num_element = _num_element;
			alphaC = new myfloat[num_element];
		}
		void free() {
			delete alphaC;
		}
		void cuda_alloc(myint _num_element) {
			num_element = _num_element;
			num_reduce = ReduceHelper::get_next_power_2_number(_num_element);
			cudaMalloc(&alphaC, num_reduce * sizeof(myfloat));
			cudaMalloc(&dev_output, ReduceHelper::get_dev_output_length(num_reduce) * sizeof(myfloat));
		}
		void cuda_free() {
			cudaFree(alphaC);
			cudaFree(dev_output);
		}
		static void cuda_memcpy(ElementFieldVariable_dt* dst, const ElementFieldVariable_dt* src, cudaMemcpyKind kind);
	
	private:
		
	};

}


#endif