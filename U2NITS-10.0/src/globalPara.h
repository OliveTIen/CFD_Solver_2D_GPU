#ifndef GLOBALPARA_H
#define GLOBALPARA_H
#include <cmath>
#include <string>


namespace Constant {
	extern const double R;//���峣��
	extern const double PI;
	extern double T0;//��ƽ���¶Ȳο�ֵ
	extern double p0;//��ƽ��ѹ���ο�ֵ
	extern double c0;//��ƽ�����ٲο�ֵ[δʹ��]
	extern double gamma;
	extern double epsilon;
}


namespace GlobalPara {

	//�ļ���
	namespace basic {
		extern int dimension;
		extern bool _continue;
		extern std::string filename;
		extern std::string meshFileType;
	}
	namespace space {
		namespace _1D {
			extern int nElement;
			extern double x1;
			extern double x2;
		}
		extern int flag_reconstruct;// �ع�������[todo]��ͨ�����췽���׻�����fun3d��δ�ҵ���ѡ��
	}
	namespace time {
		extern double CFL;
		extern double T;
		extern double residual;// �в�����
		extern int time_advance;// ʱ���ƽ���ʽ
	}
	namespace physicsModel {
		extern int equation;//"equation:1-Eluer,2-NS": 1
	}
	namespace boundaryCondition {
		namespace _2D {
			namespace inlet {
				extern bool use_ruvp;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
			namespace outlet {
				extern bool use_ruvp;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
			namespace inf {
				extern bool use_ruvp;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
		}
	}
	namespace initialCondition {
		extern int type;
	}
	namespace output {
		extern int step_per_print;
		extern int step_per_output;
		extern int step_per_output_hist;
		extern bool output_var_ruvp[4];
		extern int autosaveFileNum;
	}
	namespace inviscid_flux_method {
		extern int flux_conservation_scheme;// ��ճͨ�� �غ��ʽ ��������������� LLF Roe
		//extern int flux_construction_lhs;
		extern int flux_limiter;// ͨ��������
	}
}


#endif