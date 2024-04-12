#include "GlobalPara.h"
#include "../input/AirParameterConverter.h"

namespace GlobalPara {
	namespace constant {
		double const R = 287.06;// 气体常数
		double const PI = 3.1415926535897;
		double T0 = 288.16;// 海平面温度参考值
		double p0 = 101325.0;// 海平面压力参考值
		double c0 = 340.28;// 海平面声速参考值
		double gamma = 1.4;
		double epsilon = 1e-7;
		double Re = 1.0e8;
		double Pr = 0.73;
		double mu = 17.9e-6;// 空气动力粘性系数
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

		double CFL = 0.6;
		double CFL_steady = 0.95;
		double max_physical_time = 20;
		int time_advance = _EVO_explicit;

		// 以下数据从ContinueFile中读取，不是Toml
		double t_previous = 0;
		int istep_previous = 0;
	}
	namespace physicsModel {
		int equation = 1;//"equation:1-Eluer,2-NS": 1
	}
	namespace boundaryCondition {
		namespace _2D {
			namespace inlet {
				int input_mode = 0;
				double Ma;
				double AOA;
				double ruvp[4];
			}
			namespace outlet {
				int input_mode = 0;
				double Ma;
				double AOA;
				double ruvp[4];
			}
			namespace inf {
				int input_mode = 0;
				double Ma = 0.8;
				double AOA = 1.25;//迎角
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
		int flux_conservation_scheme = _SOL_Roe;// 黎曼求解器
		int flux_limiter = _LIM_minmod;
		int flag_reconstruct = _REC_constant;
		int flag_gradient = _GRA_leastSquare;

	}
}

