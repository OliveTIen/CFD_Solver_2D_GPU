#include "ResidualCalculator.h"
#include "../GlobalPara.h"
#include "../FVM_2D.h"

void ResidualCalculator::cal_error_isentropicVortex(double xmin, double ymin, double xmax, double ymax, 
	double chi, const double t_current, const int istep, const double cpu_time, const double* ruvp0) {
	//ruvp0：均匀流参数
	const double gamma = Constant::gamma;
	const double PI = Constant::PI;
	double xc = 0.5 * (xmin + xmax);
	double yc = 0.5 * (ymin + ymax);
	double error_norm1 = 0;//误差对体积的加权平均
	double error_norm2 = 0;//误差平方对体积的加权平均，然后开根号
	double error_max = 0;//误差的最大值
	double sumvol = 0;

	//定义lambda表达式
	//周期平移函数。若x超出[y1,y2]范围，将x平移至[y1,y2]中
	auto move_period = [](double x, double y1, double y2)->double {
		//保证y1<y2，否则swap
		if (y1 > y2) {
			double tmp = y1;
			y1 = y2;
			y2 = tmp;
		}
		else if (y1 == y2)return y1;//若[y1,y2]长度为0，则返回y1。这样避免除以0的情况
		double dy = y2 - y1;
		double xbar = x - y1;//平移至以y1为原点的坐标系
		double a = xbar / dy;
		xbar -= dy * floor(a);
		return xbar + y1;
	};

	auto elements = FVM_2D::pFVM2D->elements;
	//每个单元，计算其与精确解的误差
	for (int ie = 0; ie < elements.size(); ie++) {

		Element_T3& e = elements[ie];
		//变换坐标系。变换到初始时刻的位置。越界情况用move_period修正
		double x = move_period(e.x - t_current * ruvp0[1], xmin, xmax);
		double y = move_period(e.y - t_current * ruvp0[2], ymin, ymax);
		//以xc,yc为中心，坐标变换
		double xbar = x - xc;
		double ybar = y = yc;
		//计算rho的精确解
		double r2 = xbar * xbar + ybar * ybar;
		double dT = -(gamma - 1.) * chi * chi / (8. * gamma * PI * PI) * exp(1. - r2);
		double rho_exact = pow(ruvp0[3] + dT, 1. / (gamma - 1.));
		//rho的数值解
		double rho = e.U[0];
		//误差
		double error = abs(rho - rho_exact);
		double error2 = error * error;
		double vol = e.calArea(FVM_2D::pFVM2D);
		error_norm1 += error * vol;//1-范数误差
		error_norm2 += error2 * vol;//2-范数误差
		error_max = max(error_max, error);
		sumvol += vol;
	}
	error_norm1 /= sumvol;//  err_bar = sum(err_i * vol_i)/sum(vol_i) = sum(err_i * weight_i) 误差对体积的加权平均 
	error_norm2 = sqrt(error_norm2 / sumvol);//sqrt( sum(err_i^2 * vol_i)/sum(vol_i) ) = sqrt( sum(err_i^2 * weight_i) ) 误差平方的加权平均，然后开根号

	//写入文件
	std::ofstream outfile(GlobalStatic::exePath_withSlash + "output\\" + "error_isentropicVortex_" + GlobalPara::basic::filename + ".txt", std::ios::app);//追加模式
	outfile
		<< istep << "\t"
		<< cpu_time << "\t"
		<< error_norm1 << "\t"
		<< error_norm2 << "\t"
		<< error_max << "\n";
	outfile.close();


}

std::vector<double> ResidualCalculator::cal_residual(const std::vector<Element_T3>& elements_old, const std::vector<Element_T3>& elements, int NORM_TYPE) {
	double difference_U[4]{};// 全部初始化为0
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
			residual_U[0] += Math::max_(residual_U[0], Math::abs_(difference_U[0]));
			residual_U[1] += Math::max_(residual_U[1], Math::abs_(difference_U[1]));
			residual_U[2] += Math::max_(residual_U[2], Math::abs_(difference_U[2]));
			residual_U[3] += Math::max_(residual_U[3], Math::abs_(difference_U[3]));
		}
	}

	std::vector<double> residual_U_vector;
	residual_U_vector.resize(4);
	residual_U_vector[0] = residual_U[0];
	residual_U_vector[1] = residual_U[1];
	residual_U_vector[2] = residual_U[2];
	residual_U_vector[3] = residual_U[3];
	return residual_U_vector;
}

