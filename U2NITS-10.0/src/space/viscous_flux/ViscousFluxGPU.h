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
		C1�ǳ�����������ǰ���

		mu0			T0		S
		1.716e-5	273.15	110.4 ������վ���Ĳ���
		1.789e-5	288.16	      ���Ǻ�ƽ����� from US standard atmosphere

		C1 = mu0 / pow(T0, 1.5) * (T0 + S)
		mu = C1 * pow(T, 1.5) / (T + S)
		����վ��������õ���� C1 = 1.45793e-06������վ����Ǻ�
		UNITs��, Csthlnd = 117.0/T_inf, origianl.org��T_inf=690
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
