#include "GlobalPara.h"
#include "../input/AirParameterConverter.h"

namespace GlobalPara {
	namespace constant {
		double const R = 287.06;// ���峣��
		double const PI = 3.1415926535897;
		double T0 = 288.16;// ��ƽ���¶Ȳο�ֵ
		double p0 = 101325.0;// ��ƽ��ѹ���ο�ֵ
		double c0 = 340.28;// ��ƽ�����ٲο�ֵ
		double gamma = 1.4;
		double epsilon = 1e-7;
		double Re = 1.0e8;
		double Pr = 0.73;
		double mu = 17.9e-6;// ��������ճ��ϵ��
	}

	namespace basic {
		int dimension = 2;
		bool _continue = 1;
		std::string filename = "default_project";
		std::string meshFileType = "inp";
		bool useGPU = false;
	}
	namespace space {
		namespace _1D {
			int nElement = 300;
			double x1 = 0.0;
			double x2 = 1.0;
		}
		int flag_reconstruct = _REC_constant;
		int flag_gradient = _GRA_leastSquare;
	}
	namespace time {
		double CFL = 0.6;
		double T = 0.01;
		double residual = 1e-7;
		int time_advance = _EVO_explicit;

		// �������ݴ�ContinueFile�ж�ȡ������Toml
		double t_previous = 0;
		int istep_previous = 0;
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
		int step_per_print = 50;
		int step_per_output_field = 50;
		int step_per_output_hist = 50;
		bool output_var_ruvp[4] = { 1,1,1,1 };
		int autosaveFileNum = 3;
		int maxIteration = 5000;
	}
	namespace inviscid_flux_method {
		int flux_conservation_scheme = _SOL_Roe;// ���������
		//int flux_construction_lhs;
		int flux_limiter = _LIM_minmod;
	}
}

