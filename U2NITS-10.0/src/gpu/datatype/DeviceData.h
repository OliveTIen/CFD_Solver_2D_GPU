#ifndef DEVICE_DATA_H
#define DEVICE_DATA_H
#include "Define.h"
/*
��Щ����ֻ��Ҫ��CPU������GPU��������Ҫ������������Զ������
������CPU��һ��ָ���ʼ������Ҫ����ʱ��
*/
namespace GPU {
	//class DeviceData {
	//public:
	//	DeviceData() {}
	//	~DeviceData() {}

	//};

	/*
	DeviceArray���豸�ڴ��ϵ����顣��host��ʼ��������ʱ�Զ��ͷ��豸�ڴ档
	����Ҫʵ�ּ���Ļ�����Ҫ��д����������鷳����˽�ָ����Ϊ����
	���⣬��֪����Ա�����Ƿ���Ҫ����__device__
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
	�豸�ڴ��ϵ���������int�Ȼ�������
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
