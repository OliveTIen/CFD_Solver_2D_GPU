#include "AirParameterConverter.h"

void AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(myfloat* ruvp, myfloat Ma, myfloat AoA, myfloat T0, myfloat p0, myfloat gamma, myfloat R) {
	// ���������T0 gamma Ma p0 R
	// ����p0 T0Ϊ��ƽ�����
	// ���ݺ�ƽ��T0��gamma��Ma������T
	// ���ݺ�ƽ��p0��gamma��T��T0������p
	// ����gamma R T����c
	// ����p R T����rho
	// ����c Ma AOA����u v
	myfloat& rho = ruvp[0];
	myfloat& u = ruvp[1];
	myfloat& v = ruvp[2];
	myfloat& p = ruvp[3];
	myfloat T, c;
	myfloat const PI = 3.1415926535897;
	T = T0 / (1.0 + (gamma - 1.0) / 2.0 * Ma * Ma);
	p = p0 * pow(T / T0, gamma / (gamma - 1.0));
	c = sqrt(gamma * R * T);//��������
	rho = p / (R * T);//�����ܶ�
	u = c * Ma * cos(AoA / 180 * PI);//����x�����ٶ�
	v = c * Ma * sin(AoA / 180 * PI);//����y�����ٶ�

}
