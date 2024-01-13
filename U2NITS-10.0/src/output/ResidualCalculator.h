#ifndef RESIDUAL_CALCULATOR_H
#define RESIDUAL_CALCULATOR_H

class ResidualCalculator {
public:
	//计算等熵涡与精确解误差并写文件。ruvp0表示均匀流参数
	static void cal_error_isentropicVortex(double xmin, double ymin, double xmax, double ymax, double chi, const double t, const int istep, const double cpu_time, const double* ruvp0);
};

#endif // !RESIDUAL_CALCULATOR_H
