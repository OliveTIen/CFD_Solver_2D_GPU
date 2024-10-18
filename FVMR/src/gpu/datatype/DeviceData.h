#ifndef DEVICE_DATA_H
#define DEVICE_DATA_H
#include "DefineType.h"
#include "../../global/Constexpr.h"
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

	using DReal = DSingle<myfloat>;
	using DInt = DSingle<int>;

	class DGlobalPara {
	public:
		DArray<myfloat>* ruvp_inf;
		DArray<myfloat>* ruvp_inlet;
		DArray<myfloat>* ruvp_outlet;
		DReal* constant_R;
		DReal* constant_gamma;
		DInt* scheme;

		DGlobalPara(myfloat* inf, myfloat* inlet, myfloat* outlet, unsigned n, const myfloat* _R, const myfloat* _gamma, const int* _scheme) {
			ruvp_inf = new DArray<myfloat>(inf, n);
			ruvp_inlet = new DArray<myfloat>(inlet, n);
			ruvp_outlet = new DArray<myfloat>(outlet, n);
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


	// struct device parameter; 如果没有指针，只是传入，无需传出，就不需要申请cuda内存
	struct SDevicePara {

		/* struct type define */
		struct Constant {
			myfloat const R = 287.06;// 气体常数
			myfloat const PI = 3.1415926535897;
			myfloat T0 = 288.16;// 海平面温度参考值
			myfloat p0 = 101325.0;// 海平面压力参考值
			myfloat c0 = 340.28;// 海平面声速参考值
			myfloat gamma = 1.4;
			myfloat epsilon = 1e-7;
			myfloat Re = 1.0e8;
			myfloat Pr = 0.73;
			myfloat mu = 17.9e-6;// 空气动力粘性系数

			void initialize(myfloat _T0, myfloat _p0, myfloat _c0, myfloat _gamma, myfloat _epsilon, myfloat _Re, myfloat _Pr, myfloat _mu) {
				T0 = _T0;
				p0 = _p0;
				c0 = _c0;
				gamma = _gamma;
				epsilon = _epsilon;
				Re = _Re;
				Pr = _Pr;
				mu = _mu;
			}
		};
		struct Space {
			int flag_reconstruct = _REC_constant;
			int flag_gradient = _GRA_leastSquare;

			void initialize(int _flag_reconstruct, int _flag_gradient) {
				flag_reconstruct = _flag_reconstruct;
				flag_gradient = _flag_gradient;
			}
		};

		struct BoundaryCondition_2D_simple {
			myfloat ruvp_inf[4];
			myfloat ruvp_inlet[4];
			myfloat ruvp_outlet[4];

			void initialize(myfloat* inf, myfloat* inlet, myfloat* outlet) {
				for (int i = 0; i < 4; i++) {
					ruvp_inf[i] = inf[i];
					ruvp_inlet[i] = inlet[i];
					ruvp_outlet[i] = outlet[i];
				}
			}
		};

		struct InviscidFluxMethod {
			int flux_conservation_scheme = _SOL_Roe;// 黎曼求解器
			int flux_limiter = _LIM_minmod;

			void initialize(int _flux_conservation_scheme, int _flux_limiter) {
				flux_conservation_scheme = _flux_conservation_scheme;
				flux_limiter = _flux_limiter;
			}
		};


		/* members */
		Constant constant;
		Space space;
		BoundaryCondition_2D_simple boundaryCondition_2D;
		InviscidFluxMethod inviscidFluxMethod;

		/* functions */
		void initialize(
			myfloat _T0, myfloat _p0, myfloat _c0, myfloat _gamma, myfloat _epsilon, myfloat _Re, myfloat _Pr, myfloat _mu,
			int _flag_reconstruct, int _flag_gradient,
			myfloat* inf, myfloat* inlet, myfloat* outlet,
			int _flux_conservation_scheme, int _flux_limiter
		) {
			constant.initialize(_T0, _p0, _c0, _gamma, _epsilon, _Re, _Pr, _mu);
			space.initialize(_flag_reconstruct, _flag_gradient);
			boundaryCondition_2D.initialize(inf, inlet, outlet);
			inviscidFluxMethod.initialize(_flux_conservation_scheme, _flux_limiter);
		}
	};

}

#endif // !DEVICE_DATA_H
