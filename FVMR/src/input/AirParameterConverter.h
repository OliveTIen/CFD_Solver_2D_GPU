#ifndef AIR_PARAMETER_CONVERTER_H
#define AIR_PARAMETER_CONVERTER_H
#include <cmath>
#include "../gpu/datatype/DefineType.h"

class AirParameterConverter {
public:
	static void get_ruvp_by_Ma_AoA_T0_p0_gamma_R(myfloat* ruvp, myfloat Ma, myfloat AoA, myfloat T0, myfloat p0, myfloat gamma, myfloat R);

	static inline myfloat get_T_by_Ma_T0_gamma(myfloat Ma, myfloat T0, myfloat gamma) { return T0 / (1.0 + (gamma - 1.0) / 2.0 * Ma * Ma); }
	static inline myfloat get_p_by_p0_T_T0_gamma(myfloat p0,myfloat T,myfloat T0,myfloat gamma) { return p0 * pow(T / T0, gamma / (gamma - 1.0)); }

};
#endif // !AIR_PARAMETER_CONVERTER_H
