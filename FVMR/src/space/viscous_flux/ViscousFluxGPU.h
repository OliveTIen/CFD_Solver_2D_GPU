#ifndef VISCOUS_FLUX_GPU_H
#define VISCOUS_FLUX_GPU_H
#include "../../gpu/datatype/DefineType.h"
#include "../../gpu/Env.h"

namespace GPU {
	namespace Space {
		void edge_viscous_flux_device();

		constexpr myfloat S_Sutherland = 110.4;

		/*
		https://www.cfd-online.com/Wiki/Sutherland%27s_law
		C1是常数，可以提前算好

		mu0			T0		S
		1.716e-5	273.15	110.4 这是网站给的参数
		1.789e-5	288.16	      这是海平面参数 from US standard atmosphere

		C1 = mu0 / pow(T0, 1.5) * (T0 + S)
		mu = C1 * pow(T, 1.5) / (T + S)
		按网站参数计算得到结果 C1 = 1.45793e-06，跟网站结果吻合
		某结构程序中, Csthlnd = 117.0/T_inf, origianl.org中T_inf=690
		*/

		__host__ __device__ inline myfloat get_Sutherland_C1_host_device(myfloat S, myfloat mu0, myfloat T0) {
			return mu0 / pow(T0, 1.5) * (T0 + S);
		}

		__host__ __device__ inline myfloat get_mu_using_Sutherland_air_host_device(myfloat temperature, myfloat sutherland_C1) {
			return sutherland_C1 * pow(temperature, 1.5) / (temperature + 110.4);
		}
	}
}

#endif // !VISCOUS_FLUX_GPU_H
