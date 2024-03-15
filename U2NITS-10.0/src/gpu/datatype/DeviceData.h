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
	Device�ڴ��ϵ����顣��host��ʼ��������ʱ�Զ��ͷ��豸�ڴ档
	����Ҫʵ�ּ���Ļ�����Ҫ��д����������鷳
	���⣬��֪����Ա�����Ƿ���Ҫ����__device__
	*/
	template<typename T>
	class DArray{
	public:
		T* ptr = nullptr;
	public:
		DArray(T* v, unsigned n) {
			cudaMalloc((void**)&ptr, n * sizeof(T));
			cudaMemcpy(ptr, v, n * sizeof(T), cudaMemcpyHostToDevice);
		}
		~DArray() {
			cudaFree(ptr);
		}
	};

	class DInfPara {
	public:
		DArray<real>* ruvp_inf;
		DArray<real>* ruvp_inlet;
		DArray<real>* ruvp_outlet;

		DInfPara(real* inf, real* inlet, real* outlet, unsigned n) {
			ruvp_inf = new DArray<real>(inf, n);
			ruvp_inlet = new DArray<real>(inlet, n);
			ruvp_outlet = new DArray<real>(outlet, n);
		}
		~DInfPara() {
			delete ruvp_inf;
			delete ruvp_inlet;
			delete ruvp_outlet;
		}
	public:

	};

}

#endif // !DEVICE_DATA_H
