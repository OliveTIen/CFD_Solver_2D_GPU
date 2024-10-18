#ifndef _FIELD_SOA_H_
#define _FIELD_SOA_H_

#include "DefineType.h"
#include "../Env.h"
#include <vector>
#include "../GPUGlobalFunction.h"

namespace GPU {
	struct ReduceHelper {
	public:
		// 计算规约操作中dev_output数组长度
		static myint get_dev_output_length(myint n);

		// 获取大于等于n的最小的2的幂次方
		static myint get_next_power_2_number(myint n);
	};


	struct ElementFieldSoA {
	public:
		myint num;

		myfloat* U[4];
		myfloat* Ux[4];
		myfloat* Uy[4];
		myfloat* Flux[4];

		// 2024-05-16 添加calculateDtGPU后，ruvp在host和device均存在
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
		static void cuda_memcpy(EdgeFieldSoA* dst, const EdgeFieldSoA* src, cudaMemcpyKind kind);
	};

	/*
	用于输出Tecplot文件的场变量。
	ruvp采用传统的alloc和free，因为一定会输出到流场
	但其他变量要根据config确定是否输出，因此用动态数组

	初始化时，FieldWriter根据其OutputScheme确定vector大小
	参见FieldWriter::allocDataUsingOutputScheme
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

	struct ElementFieldVariable_dt {
		myint num_element = 0;// host模式下alphaC的长度
		myint num_reduce = 0;// 规约的数组长度，即device模式下alphaC的长度，取最近的大于num_element的2的幂
		myfloat* alphaC = nullptr;// 是临时变量，既存储alphaC，又存储dt
		myfloat* dev_output = nullptr;// 用于规约操作的临时变量。仅在device模式中申请内存，host模式禁止使用

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