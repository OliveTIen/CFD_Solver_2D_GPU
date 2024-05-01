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

	//文件类
	namespace basic {
		extern int dimension;
		extern bool _continue;
		extern std::string filename;
		extern std::string meshFileType;
		extern int useGPU;// 必须有extern，否则报错“重定义”
		extern bool isDebugMode;
	}
	//namespace space {
	//	extern int flag_reconstruct;// 重构方法。[todo]与通量构造方法易混淆，fun3d中未找到此选项
	//	extern int flag_gradient;
	//}
	namespace time {
		extern bool is_steady;
		extern bool is_explicit;

		extern double CFL;
		extern double CFL_steady;
		extern double max_physical_time;// 最大物理时间，控制终止计算
		extern int time_advance;// 时间推进方式

		// 以下数据从ContinueFile中读取，不是Toml
		extern double t_previous;// 续算时的起始时间 从readContineFile()中读取
		extern int istep_previous;// 续算时的起始步，默认为0 从readContineFile()中读取
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
		extern int flux_conservation_scheme;// 无粘通量 守恒格式 用于求解黎曼问题 LLF Roe
		//extern int flux_construction_lhs;
		extern int flux_limiter;// 通量限制器
		extern int flag_reconstruct;// 重构方法。[todo]与通量构造方法易混淆，fun3d中未找到此选项
		extern int flag_gradient;

	}
}


#endif