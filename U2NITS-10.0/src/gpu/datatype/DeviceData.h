#ifndef DEVICE_DATA_H
#define DEVICE_DATA_H
#include "Define.h"
/*
有些数据只需要从CPU拷贝到GPU，而不需要拷回来，例如远场参数
该类用CPU的一个指针初始化，需要返回时用
*/
namespace GPU {
	//class DeviceData {
	//public:
	//	DeviceData() {}
	//	~DeviceData() {}

	//};

	/*
	DeviceArray，设备内存上的数组。用host初始化，析构时自动释放设备内存。
	但是要实现计算的话，需要重写运算符，很麻烦，因此将指针设为公有
	此外，不知道成员函数是否需要加上__device__
	*/
	template<typename T>
	class DArray{
	public:
		T* ptr = nullptr;
	public:
		DArray(const T* v, unsigned n) {
			cudaMalloc((void**)&ptr, n * sizeof(T));
			cudaMemcpy(ptr, v, n * sizeof(T), cudaMemcpyHostToDevice);
		}
		~DArray() {
			cudaFree(ptr);
		}
	};

	/*
	设备内存上的数，例如int等基本类型
	*/
	template<typename T>
	class DSingle :public DArray<T> {
	public:
		DSingle(const T* v): DArray<T>(v,1) {}

	};

	using DReal = DSingle<real>;
	using DInt = DSingle<int>;

	class DGlobalPara {
	public:
		DArray<real>* ruvp_inf;
		DArray<real>* ruvp_inlet;
		DArray<real>* ruvp_outlet;
		DReal* constant_R;
		DReal* constant_gamma;
		DInt* scheme;

		DGlobalPara(real* inf, real* inlet, real* outlet, unsigned n, const real* _R, const real* _gamma, const int* _scheme) {
			ruvp_inf = new DArray<real>(inf, n);
			ruvp_inlet = new DArray<real>(inlet, n);
			ruvp_outlet = new DArray<real>(outlet, n);
			constant_R = new DReal(_R);
			constant_gamma = new DReal(_gamma);
			scheme = new DInt(_scheme);
		}
		~DGlobalPara() {
			delete ruvp_inf;
			delete ruvp_inlet;
			delete ruvp_outlet;
			delete constant_R;
			delete constant_gamma;
		}
	public:

	};

}

#endif // !DEVICE_DATA_H
