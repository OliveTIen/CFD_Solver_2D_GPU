#ifndef RESIDUAL_CALCULATOR_H
#define RESIDUAL_CALCULATOR_H

class ResidualCalculator {
public:
	//����������뾫ȷ����д�ļ���ruvp0��ʾ����������
	static void cal_error_isentropicVortex(double xmin, double ymin, double xmax, double ymax, double chi, const double t, const int istep, const double cpu_time, const double* ruvp0);
};

#endif // !RESIDUAL_CALCULATOR_H
