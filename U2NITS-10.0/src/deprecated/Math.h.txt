#ifndef M_MATH_H
#define M_MATH_H
#include <vector>
#include <iostream>
#include <string>
//#include "./include/Eigen/Core"
#include "gpu/datatype/DefineType.h"

//慎用inline，会报错LNK20219
class Math {
public:
	//---标量计算
	static inline myfloat abs_(myfloat x) { return (x > 0) ? x : (-x); }
	static inline myfloat max_(myfloat x, myfloat y) { return (x > y) ? x : y; }

	//---数组修改
	//将数组v中每个元素变为原来s倍
	static void vector_times_scalar(myfloat* v, const int vlength, const myfloat s);

	//---向量计算
	//向量点积
	static myfloat dot(std::vector<myfloat> a, std::vector<myfloat> b);
	//向量点积 仅适用于长度为2
	static inline myfloat dot2(std::vector<myfloat> a, std::vector<myfloat> b){ return a[0] * b[0] + a[1] * b[1]; }
	//向量1-范数，元素绝对值之和
	static inline myfloat vector_norm_1(std::vector<myfloat>v);
	static inline myfloat vector_norm_2(std::vector<myfloat>v);
	static inline myfloat vector_norm_infinity(std::vector<myfloat>v);


	//---时间计算
	//秒数dt转换为h:m:s
	static std::string timeFormat(int dt);
	
	//---常微分方程
	static inline myfloat RK1(myfloat Yold, myfloat h, myfloat f) { return Yold + h * f; }
	static myfloat RK3(myfloat x, myfloat y, myfloat h, myfloat (*f)(myfloat,myfloat));

	//计算当地声速
	static inline myfloat calculate_soundSpeed(myfloat gamma, myfloat p, myfloat rho) { return sqrt(gamma * p / rho); };
	static inline myfloat calculate_soundSpeed_by_temperature(myfloat gamma, myfloat R, myfloat T) { return sqrt(gamma * R * T); };
	static inline myfloat calculate_soundSpeed_by_Mach(myfloat T0, myfloat gamma, myfloat R, myfloat Ma){ 
		//常数: T0, gamma, R		//用户变量: Ma
		myfloat T = T0 / (1.0 + (gamma - 1.0) / 2.0 * Ma * Ma);//来流温度
		return calculate_soundSpeed_by_temperature(gamma, R, T);
	}

	
};

class Math_2D:public Math {
public:

	static inline myfloat get_p(const myfloat rho, const myfloat gamma, const myfloat E, const myfloat u, const myfloat v){
		return rho * (gamma - (myfloat)1.0) * (E - (u * u + v * v) * (myfloat)0.5); }
	//输入ruvp，输出rhoE
	static inline myfloat get_rhoE(const myfloat* ruvp, const myfloat gamma) {
		return ruvp[3] / (gamma - 1) + (myfloat)0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);	}//rhoE;
	static inline myfloat get_E(const myfloat* ruvp, const myfloat gamma) {
		return ruvp[3] / (gamma - 1)/ ruvp[0] + (myfloat)0.5  * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);	}//E;
	//守恒量转ruvp。
	static void U_2_ruvp(const myfloat* U, myfloat* ruvp, myfloat gamma);
	//ruvp转守恒量。
	static void ruvp_2_U(const myfloat* ruvp, myfloat* U, myfloat gamma);
	//
	static void U_2_F(const myfloat* U, myfloat* F, myfloat gamma);
	//
	static void ruvp_2_F(const myfloat* ruvp, myfloat* F, myfloat gamma);
	static void ruvp_2_Fn_lambda_2D(const myfloat* ruvp, myfloat* Fn, myfloat& lambda, myfloat nx, myfloat ny, myfloat gamma);
	static void U_2_Fn_lambda_2D(const myfloat* U, myfloat* Fn, myfloat& lambda, myfloat nx, myfloat ny, myfloat gamma);

};

#endif