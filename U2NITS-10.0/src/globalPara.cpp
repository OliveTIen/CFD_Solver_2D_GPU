#include "GlobalPara.h"
#include "SignDefine.h"
#include "input/AirParameterConverter.h"

namespace Constant {
	double const R = 287.06;//气体常数
	double const PI = 3.1415926535897;
	double T0 = 288.16;//海平面温度参考值
	double p0 = 101325.0;//海平面压力参考值
	double c0 = 340.28;//海平面声速参考值
	double gamma = 1.4;
	double epsilon = 1e-7;// 1e-7
}

namespace GlobalPara {



	namespace basic {
		int dimension = 2;
		bool _continue = 1;
		std::string filename = "default_project";
		std::string meshFileType = "inp";
	}
	namespace space {
		namespace _1D {
			int nElement = 300;
			double x1 = 0.0;
			double x2 = 1.0;
		}
		int flag_reconstruct = _REC_constant;
	}
	namespace time {
		double CFL = 0.6;
		double T = 0.01;
		double residual = 1e-7;
		int time_advance = _EVO_explicit;
	}
	namespace physicsModel {
		int equation = 1;//"equation:1-Eluer,2-NS": 1
	}
	namespace boundaryCondition {
		namespace _2D {
			namespace inlet {
				bool use_ruvp = 0;//1-使用ruvp，0-使用Ma和AOA
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
		int step_per_output = 50;
		int step_per_output_hist = 50;
		bool output_var_ruvp[4] = { 1,1,1,1 };
		int autosaveFileNum = 3;
	}
	namespace inviscid_flux_method {
		int flux_conservation_scheme = _SOL_Roe;// 黎曼求解器
		//int flux_construction_lhs;
		int flux_limiter = _LIM_minmod;
	}
}

