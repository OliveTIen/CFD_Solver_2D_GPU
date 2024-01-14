#ifndef RESIDUAL_CALCULATOR_H
#define RESIDUAL_CALCULATOR_H

#include <vector>
class Element_T3;// 声明。这样就不必包含头文件

class ResidualCalculator {
public:
	static const int NORM_1 = 1;// 1-范数，元素绝对值之和
	static const int NORM_2 = 2;// 2-范数，欧氏距离
	static const int NORM_INF = 999;// 无穷范数，元素绝对值最大值
	

public:
	// 计算等熵涡与精确解误差并写文件。ruvp0表示均匀流参数
	// ![todo]该函数需要拆分为 cal_error_xxx 和 write_xxx
	static void cal_error_isentropicVortex(double xmin, double ymin, double xmax, double ymax, double chi, const double t, const int istep, const double cpu_time, const double* ruvp0);
	static std::vector<double> cal_residual(const std::vector<Element_T3>& elements_old, const std::vector<Element_T3>& elements, int NORM_TYPE = NORM_INF);
};

#endif // !RESIDUAL_CALCULATOR_H
