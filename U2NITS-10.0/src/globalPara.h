#ifndef GLOBALPARA_H
#define GLOBALPARA_H
#include <cmath>
#include <string>


namespace GlobalPara {
	namespace constant {
		extern double const R;// ���峣��
		extern double const PI;
		extern double T0;// ��ƽ���¶Ȳο�ֵ
		extern double p0;// ��ƽ��ѹ���ο�ֵ
		extern double c0;// ��ƽ�����ٲο�ֵ
		extern double gamma;
		extern double epsilon;
		extern double Re;
		extern double Pr;
		extern double mu;// ����ճ��ϵ��
	}

	//�ļ���
	namespace basic {
		extern int dimension;
		extern bool _continue;
		extern std::string filename;
		extern std::string meshFileType;
		extern bool useGPU;// ������extern�����򱨴��ض��塱
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
		// �������ݴ�ContinueFile�ж�ȡ������Toml
		extern double t_previous;// ����ʱ����ʼʱ�� ��readContineFile()�ж�ȡ
		extern int istep_previous;// ����ʱ����ʼ�� ��readContineFile()�ж�ȡ
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
		extern int step_per_output_field;
		extern int step_per_output_hist;
		extern bool output_var_ruvp[4];
		extern int autosaveFileNum;
		extern int maxIteration;
	}
	namespace inviscid_flux_method {
		extern int flux_conservation_scheme;// ��ճͨ�� �غ��ʽ ��������������� LLF Roe
		//extern int flux_construction_lhs;
		extern int flux_limiter;// ͨ��������
	}
}


#endif