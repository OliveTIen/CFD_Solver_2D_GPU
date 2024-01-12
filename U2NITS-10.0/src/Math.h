#ifndef M_MATH_H
#define M_MATH_H
#include <vector>
#include <iostream>
#include <string>
#include "./include/Eigen/Core"

//����inline���ᱨ��LNK20219
class Math {
public:
	static double dot(std::vector<double> a, std::vector<double> b);//�������
	static inline double dot2(std::vector<double> a, std::vector<double> b){ return a[0] * b[0] + a[1] * b[1]; }
	//�������ˣ��޸�������
	static void vector_times_scalar(double* v, const int vlength, const double s);
	//����dtת��Ϊh:m:s
	static std::string timeFormat(int dt);
	//double RK1
	static inline double RK1(double Yold, double h, double f) { return Yold + h * f; }
	static double RK3(double x, double y, double h, double (*f)(double,double));
	//���㵱������
	static inline double get_c(double gamma, double p, double rho) { return sqrt(gamma * p / rho); };
	static inline double get_c_by_T(double gamma, double R, double T) { return sqrt(gamma * R * T); };
	static inline double get_c_by_Ma(double T0, double gamma, double R, double Ma){ 
		//����: T0, gamma, R		//�û�����: Ma
		double T = T0 / (1.0 + (gamma - 1.0) / 2.0 * Ma * Ma);//�����¶�
		return sqrt(gamma * R * T);//��������
	}
};

class Math_2D:public Math {
public:

	static inline double get_p(const double rho, const double gamma, const double E, const double u, const double v){
		return rho * (gamma - 1) * (E - (u * u + v * v) * 0.5); }
	//����ruvp�����rhoE
	static inline double get_rhoE(const double* ruvp, const double gamma) {
		return ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);	}//rhoE;
	static inline double get_E(const double* ruvp, const double gamma) {
		return ruvp[3] / (gamma - 1)/ ruvp[0] + 0.5  * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);	}//E;
	//�غ���תruvp��
	static void U_2_ruvp(const double* U, double* ruvp, double gamma);
	//ruvpת�غ�����
	static void ruvp_2_U(const double* ruvp, double* U, double gamma);
	//
	static void U_2_F(const double* U, double* F, double gamma);
	static void U_2_F(const Eigen::Vector4d& U, Eigen::Vector4d& F, double gamma);
	static void U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, double& lambda, double gamma);
	//
	static void ruvp_2_F(const double* ruvp, double* F, double gamma);
	static void ruvp_2_Fn_lambda_2D(const double* ruvp, double* Fn, double& lambda, double nx, double ny, double gamma);
	static void U_2_Fn_lambda_2D(const double* U, double* Fn, double& lambda, double nx, double ny, double gamma);

};

#endif