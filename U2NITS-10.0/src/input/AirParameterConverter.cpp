#include "AirParameterConverter.h"

void AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(double* ruvp, double Ma, double AoA, double T0, double p0, double gamma, double R) {
	// 输入参数：T0 gamma Ma p0 R
	// 其中p0 T0为海平面参数
	// 根据海平面T0、gamma、Ma，计算T
	// 根据海平面p0、gamma、T、T0，计算p
	// 根据gamma R T计算c
	// 根据p R T计算rho
	// 根据c Ma AOA计算u v
	double& rho = ruvp[0];
	double& u = ruvp[1];
	double& v = ruvp[2];
	double& p = ruvp[3];
	double T, c;
	double const PI = 3.1415926535897;
	T = T0 / (1.0 + (gamma - 1.0) / 2.0 * Ma * Ma);
	p = p0 * pow(T / T0, gamma / (gamma - 1.0));
	c = sqrt(gamma * R * T);//来流声速
	rho = p / (R * T);//来流密度
	u = c * Ma * cos(AoA / 180 * PI);//来流x方向速度
	v = c * Ma * sin(AoA / 180 * PI);//来流y方向速度

}
