#ifndef GLOBALPARA_H
#define GLOBALPARA_H
#include <cmath>
#include <string>


namespace Constant {
	extern const double R;//气体常数
	extern const double PI;
	extern double T0;//海平面温度参考值
	extern double p0;//海平面压力参考值
	extern double c0;//海平面声速参考值
	extern double gamma;
}


namespace GlobalPara {

	//文件类
	namespace basic {
		extern int dimension;
		extern bool _continue;
		extern std::string filename;
	}
	namespace space {
		namespace _1D {
			extern int nElement;
			extern double x1;
			extern double x2;
		}
		extern int flag_reconstruct;// 重构
	}
	namespace time {
		extern double CFL;
		extern double T;
		extern double residual;// 残差限制
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
			//若不给定ruvp，则根据Ma和AOA计算ruvp
			void ini_ruvp_by_Ma_AOA();
			//根据海平面参数，计算ruvp
			void get_ruvp_sealevel(const double Ma, const double AOA, double* ruvp);
			inline double get_T_sealevel(const double Ma){ return Constant::T0 / (1.0 + (Constant::gamma - 1.0) / 2.0 * Ma * Ma); }
			inline double get_p_sealevel(const double T){ return Constant::p0 * pow(T / Constant::T0, Constant::gamma / (Constant::gamma - 1.0)); }
		}
	}
	namespace initialCondition {
		extern int type;
	}
	namespace output {
		extern int step_per_output;
		extern bool output_var_ruvp[4];
		extern int autosaveFileNum;
	}
}


#endif