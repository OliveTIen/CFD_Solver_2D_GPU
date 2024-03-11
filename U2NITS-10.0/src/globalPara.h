#ifndef GLOBALPARA_H
#define GLOBALPARA_H
#include <cmath>
#include <string>


namespace GlobalPara {
	namespace constant {
		extern double const R;// 气体常数
		extern double const PI;
		extern double T0;// 海平面温度参考值
		extern double p0;// 海平面压力参考值
		extern double c0;// 海平面声速参考值
		extern double gamma;
		extern double epsilon;
		extern double Re;
		extern double Pr;
		extern double mu;// 动力粘性系数
	}

	//文件类
	namespace basic {
		extern int dimension;
		extern bool _continue;
		extern std::string filename;
		extern std::string meshFileType;
		extern bool useGPU;// 必须有extern，否则报错“重定义”
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
		// 以下数据从ContinueFile中读取，不是Toml
		extern double t_previous;// 续算时的起始时间 从readContineFile()中读取
		extern int istep_previous;// 续算时的起始步 从readContineFile()中读取
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
		extern int flux_conservation_scheme;// 无粘通量 守恒格式 用于求解黎曼问题 LLF Roe
		//extern int flux_construction_lhs;
		extern int flux_limiter;// 通量限制器
	}
}


#endif