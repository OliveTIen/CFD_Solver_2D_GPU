#ifndef RESIDUAL_CALCULATOR_H
#define RESIDUAL_CALCULATOR_H

#include <vector>
class Element_2D;// �����������Ͳ��ذ���ͷ�ļ�

class ResidualCalculator {
public:
	static const int NORM_1 = 1;// 1-������Ԫ�ؾ���ֵ֮��
	static const int NORM_2 = 2;// 2-������ŷ�Ͼ���
	static const int NORM_INF = 999;// �������Ԫ�ؾ���ֵ���ֵ
	

public:
	// ����������뾫ȷ����д�ļ���ruvp0��ʾ����������
	// ![todo]�ú�����Ҫ���Ϊ cal_error_xxx �� write_xxx
	static void cal_error_isentropicVortex(double xmin, double ymin, double xmax, double ymax, double chi, const double t, const int istep, const double cpu_time, const double* ruvp0);
	static void cal_residual(const std::vector<Element_2D>& elements_old, const std::vector<Element_2D>& elements, int NORM_TYPE, double* residual);
};

#endif // !RESIDUAL_CALCULATOR_H
