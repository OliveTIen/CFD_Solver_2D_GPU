#ifndef M_MATH_H
#define M_MATH_H
#include <vector>
#include <iostream>
#include <string>
#include "./include/Eigen/Core"

//慎用inline，会报错LNK20219
class Math {
public:
	//---标量计算
	static inline double abs_(double x) { return (x > 0) ? x : (-x); }
	static inline double max_(double x, double y) { return (x > y) ? x : y; }

	//---数组修改
	//将数组v中每个元素变为原来s倍
	static void vector_times_scalar(double* v, const int vlength, const double s);

	//---向量计算
	//向量点积
	static double dot(std::vector<double> a, std::vector<double> b);
	//向量点积 仅适用于长度为2
	static inline double dot2(std::vector<double> a, std::vector<double> b){ return a[0] * b[0] + a[1] * b[1]; }
	//向量1-范数，元素绝对值之和
	static inline double vector_norm_1(std::vector<double>v);
	static inline double vector_norm_2(std::vector<double>v);
	static inline double vector_norm_infinity(std::vector<double>v);


	//---时间计算
	//秒数dt转换为h:m:s
	static std::string timeFormat(int dt);
	
	//---常微分方程
	static inline double RK1(double Yold, double h, double f) { return Yold + h * f; }
	static double RK3(double x, double y, double h, double (*f)(double,double));

	//计算当地声速
	static inline double calculate_soundSpeed(double gamma, double p, double rho) { return sqrt(gamma * p / rho); };
	static inline double calculate_soundSpeed_by_temperature(double gamma, double R, double T) { return sqrt(gamma * R * T); };
	static inline double calculate_soundSpeed_by_Mach(double T0, double gamma, double R, double Ma){ 
		//常数: T0, gamma, R		//用户变量: Ma
		double T = T0 / (1.0 + (gamma - 1.0) / 2.0 * Ma * Ma);//来流温度
		return calculate_soundSpeed_by_temperature(gamma, R, T);
	}

	
};

class Math_2D:public Math {
public:

	static inline double get_p(const double rho, const double gamma, const double E, const double u, const double v){
		return rho * (gamma - 1) * (E - (u * u + v * v) * 0.5); }
	//输入ruvp，输出rhoE
	static inline double get_rhoE(const double* ruvp, const double gamma) {
		return ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);	}//rhoE;
	static inline double get_E(const double* ruvp, const double gamma) {
		return ruvp[3] / (gamma - 1)/ ruvp[0] + 0.5  * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);	}//E;
	//守恒量转ruvp。
	static void U_2_ruvp(const double* U, double* ruvp, double gamma);
	//ruvp转守恒量。
	static void ruvp_2_U(const double* ruvp, double* U, double gamma);
	//
	static void U_2_F(const double* U, double* F, double gamma);
	//
	static void ruvp_2_F(const double* ruvp, double* F, double gamma);
	static void ruvp_2_Fn_lambda_2D(const double* ruvp, double* Fn, double& lambda, double nx, double ny, double gamma);
	static void U_2_Fn_lambda_2D(const double* U, double* Fn, double& lambda, double nx, double ny, double gamma);

};

#endif