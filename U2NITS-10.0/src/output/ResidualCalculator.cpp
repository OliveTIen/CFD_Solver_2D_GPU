#include "ResidualCalculator.h"
#include "../GlobalPara.h"
#include "../FVM_2D.h"

void ResidualCalculator::cal_error_isentropicVortex(double xmin, double ymin, double xmax, double ymax, double chi, const double t_current, const int istep, const double cpu_time, const double* ruvp0) {
	//ruvp0������������
	const double gamma = Constant::gamma;
	const double PI = Constant::PI;
	double xc = 0.5 * (xmin + xmax);
	double yc = 0.5 * (ymin + ymax);
	double error_norm1 = 0;//��������ļ�Ȩƽ��
	double error_norm2 = 0;//���ƽ��������ļ�Ȩƽ����Ȼ�󿪸���
	double error_max = 0;//�������ֵ
	double sumvol = 0;

	//����lambda���ʽ
	//����ƽ�ƺ�������x����[y1,y2]��Χ����xƽ����[y1,y2]��
	auto move_period = [](double x, double y1, double y2)->double {
		//��֤y1<y2������swap
		if (y1 > y2) {
			double tmp = y1;
			y1 = y2;
			y2 = tmp;
		}
		else if (y1 == y2)return y1;//��[y1,y2]����Ϊ0���򷵻�y1�������������0�����
		double dy = y2 - y1;
		double xbar = x - y1;//ƽ������y1Ϊԭ�������ϵ
		double a = xbar / dy;
		xbar -= dy * floor(a);
		return xbar + y1;
	};

	auto elements = FVM_2D::pFVM2D->elements;
	//ÿ����Ԫ���������뾫ȷ������
	for (int ie = 0; ie < elements.size(); ie++) {

		Element_T3& e = elements[ie];
		//�任����ϵ���任����ʼʱ�̵�λ�á�Խ�������move_period����
		double x = move_period(e.x - t_current * ruvp0[1], xmin, xmax);
		double y = move_period(e.y - t_current * ruvp0[2], ymin, ymax);
		//��xc,ycΪ���ģ�����任
		double xbar = x - xc;
		double ybar = y = yc;
		//����rho�ľ�ȷ��
		double r2 = xbar * xbar + ybar * ybar;
		double dT = -(gamma - 1.) * chi * chi / (8. * gamma * PI * PI) * exp(1. - r2);
		double rho_exact = pow(ruvp0[3] + dT, 1. / (gamma - 1.));
		//rho����ֵ��
		double rho = e.U[0];
		//���
		double error = abs(rho - rho_exact);
		double error2 = error * error;
		double vol = e.calArea(FVM_2D::pFVM2D);
		error_norm1 += error * vol;//1-�������
		error_norm2 += error2 * vol;//2-�������
		error_max = max(error_max, error);
		sumvol += vol;
	}
	error_norm1 /= sumvol;//  err_bar = sum(err_i * vol_i)/sum(vol_i) = sum(err_i * weight_i) ��������ļ�Ȩƽ�� 
	error_norm2 = sqrt(error_norm2 / sumvol);//sqrt( sum(err_i^2 * vol_i)/sum(vol_i) ) = sqrt( sum(err_i^2 * weight_i) ) ���ƽ���ļ�Ȩƽ����Ȼ�󿪸���

	//д���ļ�
	std::ofstream outfile(GlobalStatic::exePath_withSlash + "output\\" + "error_isentropicVortex_" + GlobalStatic::filename + ".txt", std::ios::app);//׷��ģʽ
	outfile
		<< istep << "\t"
		<< cpu_time << "\t"
		<< error_norm1 << "\t"
		<< error_norm2 << "\t"
		<< error_max << "\n";
	outfile.close();


}
