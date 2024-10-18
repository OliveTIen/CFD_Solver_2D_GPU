#include "Constexpr.h"

#ifndef GLOBALPARA_H
#define GLOBALPARA_H
#include <cmath>
#include <string>
#include "../gpu/datatype/DefineType.h"


namespace GlobalPara {
	namespace constant {
		extern myfloat R;
		extern myfloat T0;
		extern myfloat p0;
		extern myfloat c0;
		extern myfloat gamma;
		extern myfloat epsilon;
		extern myfloat Re;
		extern myfloat Pr;
		extern myfloat mu0;

		extern myfloat referenceArea;
		extern myfloat mesh_scale_factor;
	}

	namespace render {
		extern myfloat range_1;
		extern myfloat range_2;
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
	namespace time {
		extern bool is_steady;
		extern bool is_explicit;

		extern myfloat CFL;
		extern myfloat CFL_steady;
		extern myfloat max_physical_time;// 最大物理时间，控制终止计算
		extern int time_advance;// 时间推进方式

		// 以下数据从ContinueFile中读取，不是Toml
		extern myfloat t_previous;// 续算时的起始时间 从readContineFile()中读取
		extern int istep_previous;// 续算时的起始步，默认为0 从readContineFile()中读取
	}
	namespace physicsModel {
		extern int equation;//"equation:1-Eluer,2-NS": 1
	}
	namespace boundaryCondition {
		namespace _2D {
			namespace inlet {
				extern int input_mode;
				extern myfloat Ma;
				extern myfloat AOA;
				extern myfloat ruvp[4];
			}
			namespace outlet {
				extern int input_mode;
				extern myfloat Ma;
				extern myfloat AOA;
				extern myfloat ruvp[4];
			}
			namespace inf {
				extern int input_mode;
				extern myfloat Ma;
				extern myfloat AOA;
				extern myfloat ruvp[4];
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
		extern int start_output_field;
		extern int autosaveFileNum;
		extern int maxIteration;
		extern myfloat tolerace_residual;
	}
	namespace inviscid_flux_method {
		extern int flux_conservation_scheme;// 无粘通量 守恒格式 用于求解黎曼问题 LLF Roe
		//extern int flux_construction_lhs;
		extern int flux_limiter;// 通量限制器
		extern int flag_reconstruct;// 重构方法
		extern int flag_gradient;

	}
}


#endif