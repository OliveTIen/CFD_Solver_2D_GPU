#include "ResidualCalculator.h"
#include "../GlobalPara.h"
#include "../FVM_2D.h"
#include "../global/FilePathManager.h"

void ResidualCalculator::cal_error_isentropicVortex(double xmin, double ymin, double xmax, double ymax, 
	double chi, const double t_current, const int istep, const double cpu_time, const double* ruvp0) {
	//ruvp0������������
	const double gamma = GlobalPara::constant::gamma;
	const double PI = GlobalPara::constant::PI;
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

	auto elements = FVM_2D::getInstance()->elements;
	//ÿ����Ԫ���������뾫ȷ������
	for (int ie = 0; ie < elements.size(); ie++) {

		Element_2D& e = elements[ie];
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
		double vol = e.calArea(FVM_2D::getInstance());
		error_norm1 += error * vol;//1-�������
		error_norm2 += error2 * vol;//2-�������
		error_max = max(error_max, error);
		sumvol += vol;
	}
	error_norm1 /= sumvol;//  err_bar = sum(err_i * vol_i)/sum(vol_i) = sum(err_i * weight_i) ��������ļ�Ȩƽ�� 
	error_norm2 = sqrt(error_norm2 / sumvol);//sqrt( sum(err_i^2 * vol_i)/sum(vol_i) ) = sqrt( sum(err_i^2 * weight_i) ) ���ƽ���ļ�Ȩƽ����Ȼ�󿪸���

	//д���ļ�
	std::ofstream outfile(FilePathManager::getInstance()->getOutputDirectory() + "error_isentropicVortex_" + GlobalPara::basic::filename + ".txt", std::ios::app);//׷��ģʽ
	outfile
		<< istep << "\t"
		<< cpu_time << "\t"
		<< error_norm1 << "\t"
		<< error_norm2 << "\t"
		<< error_max << "\n";
	outfile.close();


}

void ResidualCalculator::cal_residual(const std::vector<Element_2D>& elements_old, const std::vector<Element_2D>& elements, int NORM_TYPE, double* residual_to_be_changed) {
	double difference_U[4]{};// ȫ����ʼ��Ϊ0
	double residual_U[4]{};

	switch (NORM_TYPE) {
	case NORM_1:
		for (int ie = 0; ie < elements.size(); ie++) {
			difference_U[0] = elements[ie].U[0] - elements_old[ie].U[0];
			difference_U[1] = elements[ie].U[1] - elements_old[ie].U[1];
			difference_U[2] = elements[ie].U[2] - elements_old[ie].U[2];
			difference_U[3] = elements[ie].U[3] - elements_old[ie].U[3];
			residual_U[0] += Math::abs_(difference_U[0]);
			residual_U[1] += Math::abs_(difference_U[1]);
			residual_U[2] += Math::abs_(difference_U[2]);
			residual_U[3] += Math::abs_(difference_U[3]);
		}
		break;
	case NORM_2:
		for (int ie = 0; ie < elements.size(); ie++) {
			difference_U[0] = elements[ie].U[0] - elements_old[ie].U[0];
			difference_U[1] = elements[ie].U[1] - elements_old[ie].U[1];
			difference_U[2] = elements[ie].U[2] - elements_old[ie].U[2];
			difference_U[3] = elements[ie].U[3] - elements_old[ie].U[3];
			residual_U[0] += difference_U[0] * difference_U[0];
			residual_U[1] += difference_U[1] * difference_U[1];
			residual_U[2] += difference_U[2] * difference_U[2];
			residual_U[3] += difference_U[3] * difference_U[3];
		}
		residual_U[0] = sqrt(residual_U[0]);
		residual_U[1] = sqrt(residual_U[1]);
		residual_U[2] = sqrt(residual_U[2]);
		residual_U[3] = sqrt(residual_U[3]);
		break;
	default:// NORM_INF
		for (int ie = 0; ie < elements.size(); ie++) {
			difference_U[0] = elements[ie].U[0] - elements_old[ie].U[0];
			difference_U[1] = elements[ie].U[1] - elements_old[ie].U[1];
			difference_U[2] = elements[ie].U[2] - elements_old[ie].U[2];
			difference_U[3] = elements[ie].U[3] - elements_old[ie].U[3];
			residual_U[0] = Math::max_(residual_U[0], Math::abs_(difference_U[0]));
			residual_U[1] = Math::max_(residual_U[1], Math::abs_(difference_U[1]));
			residual_U[2] = Math::max_(residual_U[2], Math::abs_(difference_U[2]));
			residual_U[3] = Math::max_(residual_U[3], Math::abs_(difference_U[3]));
		}
	}
	// ���
	residual_to_be_changed[0] = residual_U[0];
	residual_to_be_changed[1] = residual_U[1];
	residual_to_be_changed[2] = residual_U[2];
	residual_to_be_changed[3] = residual_U[3];
}

void ResidualCalculator::cal_residual_GPU(REAL* elements_U_old[4], GPU::FieldSoA elements, int NORM_TYPE, double* residual) {
	// ����в�
	double difference_U[4]{};// ȫ����ʼ��Ϊ0
	double residual_U[4]{};
	const int numElements = elements.num;

	switch (NORM_TYPE) {
	case NORM_1:
		for (int ie = 0; ie < numElements; ie++) {
			difference_U[0] = elements.U[0][ie] - elements_U_old[0][ie];
			difference_U[1] = elements.U[1][ie] - elements_U_old[1][ie];
			difference_U[2] = elements.U[2][ie] - elements_U_old[2][ie];
			difference_U[3] = elements.U[3][ie] - elements_U_old[3][ie];
			residual_U[0] += Math::abs_(difference_U[0]);
			residual_U[1] += Math::abs_(difference_U[1]);
			residual_U[2] += Math::abs_(difference_U[2]);
			residual_U[3] += Math::abs_(difference_U[3]);
		}
		break;
	case NORM_2:
		for (int ie = 0; ie < numElements; ie++) {
			difference_U[0] = elements.U[0][ie] - elements_U_old[0][ie];
			difference_U[1] = elements.U[1][ie] - elements_U_old[1][ie];
			difference_U[2] = elements.U[2][ie] - elements_U_old[2][ie];
			difference_U[3] = elements.U[3][ie] - elements_U_old[3][ie];
			residual_U[0] += difference_U[0] * difference_U[0];
			residual_U[1] += difference_U[1] * difference_U[1];
			residual_U[2] += difference_U[2] * difference_U[2];
			residual_U[3] += difference_U[3] * difference_U[3];
		}
		residual_U[0] = sqrt(residual_U[0]);
		residual_U[1] = sqrt(residual_U[1]);
		residual_U[2] = sqrt(residual_U[2]);
		residual_U[3] = sqrt(residual_U[3]);
		break;
	default:// NORM_INF
		for (int ie = 0; ie < numElements; ie++) {
			difference_U[0] = elements.U[0][ie] - elements_U_old[0][ie];
			difference_U[1] = elements.U[1][ie] - elements_U_old[1][ie];
			difference_U[2] = elements.U[2][ie] - elements_U_old[2][ie];
			difference_U[3] = elements.U[3][ie] - elements_U_old[3][ie];
			residual_U[0] = Math::max_(residual_U[0], Math::abs_(difference_U[0]));
			residual_U[1] = Math::max_(residual_U[1], Math::abs_(difference_U[1]));
			residual_U[2] = Math::max_(residual_U[2], Math::abs_(difference_U[2]));
			residual_U[3] = Math::max_(residual_U[3], Math::abs_(difference_U[3]));
		}
	}
	// ���
	residual[0] = residual_U[0];
	residual[1] = residual_U[1];
	residual[2] = residual_U[2];
	residual[3] = residual_U[3];

}

