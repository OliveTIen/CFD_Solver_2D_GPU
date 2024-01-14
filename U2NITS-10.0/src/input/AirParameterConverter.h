#ifndef AIR_PARAMETER_CONVERTER_H
#define AIR_PARAMETER_CONVERTER_H
#include <cmath>

class AirParameterConverter {
public:
	static void get_ruvp_by_Ma_AoA_T0_p0_gamma_R(double* ruvp, double Ma, double AoA, double T0, double p0, double gamma, double R);

	static inline double get_T_by_Ma_T0_gamma(double Ma, double T0, double gamma) { return T0 / (1.0 + (gamma - 1.0) / 2.0 * Ma * Ma); }
	static inline double get_p_by_p0_T_T0_gamma(double p0,double T,double T0,double gamma) { return p0 * pow(T / T0, gamma / (gamma - 1.0)); }

};
#endif // !AIR_PARAMETER_CONVERTER_H
