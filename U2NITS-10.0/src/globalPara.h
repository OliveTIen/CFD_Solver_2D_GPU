#ifndef GLOBALPARA_H
#define GLOBALPARA_H
#include <cmath>
#include <string>


namespace Constant {
	extern const double R;//气体常数
	extern const double PI;
	extern double T0;//海平面温度参考值
	extern double p0;//海平面压力参考值
	extern double c0;//海平面声速参考值[未使用]
	extern double gamma;
	extern double epsilon;
}


namespace GlobalPara {

	//文件类
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
		extern int flag_reconstruct;// 重构方法。[todo]与通量构造方法易混淆，fun3d中未找到此选项
	}
	namespace time {
		extern double CFL;
		extern double T;
		extern double residual;// 残差限制
		extern int time_advance;// 时间推进方式
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
		extern int flux_conservation_scheme;// 无粘通量 守恒格式 用于求解黎曼问题 LLF Roe
		//extern int flux_construction_lhs;
		extern int flux_limiter;// 通量限制器
	}
}


#endif