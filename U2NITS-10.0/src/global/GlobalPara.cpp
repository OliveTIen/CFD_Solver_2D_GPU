#include "GlobalPara.h"
#include "../input/AirParameterConverter.h"

namespace GlobalPara {
	namespace constant {
		myfloat R = 287.06;                   // ���峣��
		myfloat T0 = 288.16;                  // ��ƽ���¶Ȳο�ֵ
		myfloat p0 = 101325.0;                // ��ƽ��ѹ���ο�ֵ
		myfloat c0 = 340.28;                  // ��ƽ�����ٲο�ֵ
		myfloat gamma = 1.4;
		myfloat epsilon = 1e-7;
		myfloat Re = 1.0e8;
		myfloat Pr = 0.73;
		myfloat mu0 = 17.9e-6;                 // ��������ճ��ϵ���ο�ֵ

		myfloat referenceArea = 1.0;          // (��ά����)�ο��������(��ά����)�ο��ҳ���������ǰ��
		myfloat mesh_scale_factor = 1.0;            // �����������ӡ�ֻ�ڶ�ȡ�����ļ�ʱ��Ч���ο����Ҳ����Ӱ��
	}

	namespace basic {
		int dimension = 2;
		bool _continue = 1;
		std::string filename = "default_project";
		std::string meshFileType = "inp";
		int useGPU = 0;
		bool isDebugMode = false;
	}
	namespace time {
		bool is_steady = 0;
		bool is_explicit = 1;

		myfloat CFL = 0.6;
		myfloat CFL_steady = 0.95;
		myfloat max_physical_time = 20;
		int time_advance = _EVO_explicit;

		// �������ݴ�ContinueFile�ж�ȡ������Toml
		myfloat t_previous = 0;
		int istep_previous = 0;
	}
	namespace physicsModel {
		int equation = 1;//"equation:1-Eluer,2-NS": 1
	}
	namespace boundaryCondition {
		namespace _2D {
			namespace inlet {
				int input_mode = 0;
				myfloat Ma;
				myfloat AOA;
				myfloat ruvp[4];
			}
			namespace outlet {
				int input_mode = 0;
				myfloat Ma;
				myfloat AOA;
				myfloat ruvp[4];
			}
			namespace inf {
				int input_mode = 0;
				myfloat Ma = 0.8;
				myfloat AOA = 1.25;//ӭ��
				myfloat ruvp[4];
			}
		}
	}
	namespace initialCondition {
		//int type = 1;
	}
	namespace output {
		int step_per_print = 50;
		int step_per_output_field = 50;
		int step_per_output_hist = 50;
		int start_output_field = 0;// ��ʼ����Ĳ���
		int autosaveFileNum = 3;
		int maxIteration = 5000;
		myfloat tolerace_residual = 1.0e-7;// ���ڴ˲в���Ϊ�ﵽ��̬
	}
	namespace inviscid_flux_method {
		int flux_conservation_scheme = _SOL_Roe;// ���������
		int flux_limiter = _LIM_minmod;
		int flag_reconstruct = _REC_constant;
		int flag_gradient = _GRA_leastSquare;

	}
}

