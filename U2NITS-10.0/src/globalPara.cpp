#include "globalPara.h"
#include "SignDefine.h"

namespace Constant {
	double const R = 287.06;//���峣��
	double const PI = 3.1415926535897;
	double T0 = 288.16;//��ƽ���¶Ȳο�ֵ
	double p0 = 101325.0;//��ƽ��ѹ���ο�ֵ
	double c0 = 340.28;//��ƽ�����ٲο�ֵ
	double gamma = 1.4;
}

namespace GlobalPara {
	namespace flag {
		//int equation = 1;
		int reconstruct = _REC_constant;
	}



	namespace basic {
		int dimension = 2;
		bool _continue = 1;
	}
	namespace space {
		namespace _1D {
			int nElement = 300;
			double x1 = 0.0;
			double x2 = 1.0;
		}
	}
	namespace time {
		double CFL = 0.6;
		double T = 0.01;
	}
	namespace physicsModel {
		int equation = 1;//"equation:1-Eluer,2-NS": 1
	}
	namespace boundaryCondition {
		namespace _2D {
			namespace inlet {
				bool use_ruvp = 0;//1-ʹ��ruvp��0-ʹ��Ma��AOA
				double Ma;
				double AOA;
				double ruvp[4];
			}
			namespace outlet {
				bool use_ruvp = 0;
				double Ma;
				double AOA;
				double ruvp[4];
			}
			namespace inf {
				bool use_ruvp = 0;
				double Ma = 0.8;
				double AOA = 1.25;//ӭ��
				double ruvp[4];
			}
		}
	}
	namespace initialCondition {
		int type = 1;
	}
	namespace output {
		int step_per_output = 50;
		bool output_var_ruvp[4] = { 1,1,1,1 };
	}
}





void GlobalPara::boundaryCondition::_2D::ini_ruvp_by_Ma_AOA() {
	//��������ruvp�������Ma��AOA����ruvp
	using namespace GlobalPara::boundaryCondition::_2D;
	if (!inf::use_ruvp) 	get_ruvp_sealevel(inf::Ma, inf::AOA, inf::ruvp);
	if (!inlet::use_ruvp)	get_ruvp_sealevel(inlet::Ma, inlet::AOA, inlet::ruvp);
	if (!outlet::use_ruvp)	get_ruvp_sealevel(outlet::Ma, outlet::AOA, outlet::ruvp);
}

void GlobalPara::boundaryCondition::_2D::get_ruvp_sealevel(const double Ma, const double AOA, double* ruvp) {
	double& rho = ruvp[0];
	double& u = ruvp[1];
	double& v = ruvp[2];
	double& p = ruvp[3];
	double T, c;
	T = get_T_sealevel(Ma);//�����¶�
	p = get_p_sealevel(T);//����ѹ��
	c = sqrt(Constant::gamma * Constant::R * T);//��������
	rho = p / (Constant::R * T);//�����ܶ�
	u = c * Ma * cos(AOA / 180 * Constant::PI);//����x�����ٶ�
	v = c * Ma * sin(AOA / 180 * Constant::PI);//����y�����ٶ�
}

//double GlobalPara::boundaryCondition::_2D::get_T_sealevel(const double Ma) {
//	return Constant::T0 / (1.0 + (Constant::gamma - 1.0) / 2.0 * Ma * Ma); 
//}

