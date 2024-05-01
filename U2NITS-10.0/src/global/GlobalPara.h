#include "Constexpr.h"

#ifndef GLOBALPARA_H
#define GLOBALPARA_H
#include <cmath>
#include <string>



namespace GlobalPara {
	namespace constant {
		extern double R;
		extern double T0;
		extern double p0;
		extern double c0;
		extern double gamma;
		extern double epsilon;
		extern double Re;
		extern double Pr;
		extern double mu;

		extern double referenceArea;
	}

	//�ļ���
	namespace basic {
		extern int dimension;
		extern bool _continue;
		extern std::string filename;
		extern std::string meshFileType;
		extern int useGPU;// ������extern�����򱨴��ض��塱
		extern bool isDebugMode;
	}
	//namespace space {
	//	extern int flag_reconstruct;// �ع�������[todo]��ͨ�����췽���׻�����fun3d��δ�ҵ���ѡ��
	//	extern int flag_gradient;
	//}
	namespace time {
		extern bool is_steady;
		extern bool is_explicit;

		extern double CFL;
		extern double CFL_steady;
		extern double max_physical_time;// �������ʱ�䣬������ֹ����
		extern int time_advance;// ʱ���ƽ���ʽ

		// �������ݴ�ContinueFile�ж�ȡ������Toml
		extern double t_previous;// ����ʱ����ʼʱ�� ��readContineFile()�ж�ȡ
		extern int istep_previous;// ����ʱ����ʼ����Ĭ��Ϊ0 ��readContineFile()�ж�ȡ
	}
	namespace physicsModel {
		extern int equation;//"equation:1-Eluer,2-NS": 1
	}
	namespace boundaryCondition {
		namespace _2D {
			namespace inlet {
				extern int input_mode;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
			namespace outlet {
				extern int input_mode;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
			namespace inf {
				extern int input_mode;
				extern double Ma;
				extern double AOA;
				extern double ruvp[4];
			}
		}
	}
	namespace initialCondition {
		//extern int type;
	}
	namespace output {
		extern int step_per_print;
		extern int step_per_output_field;
		extern int step_per_output_hist;
		extern int autosaveFileNum;
		extern int maxIteration;
		extern double tolerace_residual;
	}
	namespace inviscid_flux_method {
		extern int flux_conservation_scheme;// ��ճͨ�� �غ��ʽ ��������������� LLF Roe
		//extern int flux_construction_lhs;
		extern int flux_limiter;// ͨ��������
		extern int flag_reconstruct;// �ع�������[todo]��ͨ�����췽���׻�����fun3d��δ�ҵ���ѡ��
		extern int flag_gradient;

	}
}


#endif