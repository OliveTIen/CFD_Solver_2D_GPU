#include "GlobalPara.h"
#include "../input/AirParameterConverter.h"

namespace GlobalPara {
	namespace constant {
		myfloat R = 287.06;                   // 气体常数
		myfloat T0 = 288.16;                  // 海平面温度参考值
		myfloat p0 = 101325.0;                // 海平面压力参考值
		myfloat c0 = 340.28;                  // 海平面声速参考值
		myfloat gamma = 1.4;
		myfloat epsilon = 1e-7;
		myfloat Re = 1.0e8;
		myfloat Pr = 0.73;
		myfloat mu0 = 17.9e-6;                 // 空气动力粘性系数参考值

		myfloat referenceArea = 1.0;          // (三维算例)参考面积，或(二维算例)参考弦长。填缩放前的
		myfloat mesh_scale_factor = 1.0;            // 网格缩放因子。只在读取网格文件时生效。参考面积也受其影响
	}

	namespace render {
		myfloat range_1 = 0;
		myfloat range_2 = 1;
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

		// 以下数据从ContinueFile中读取，不是Toml
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
				myfloat AOA = 1.25;//迎角
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
		int start_output_field = 0;// 开始输出的步数
		int autosaveFileNum = 3;
		int maxIteration = 5000;
		myfloat tolerace_residual = 1.0e-7;// 低于此残差认为达到稳态
	}
	namespace inviscid_flux_method {
		int flux_conservation_scheme = _SOL_Roe;// 黎曼求解器
		int flux_limiter = _LIM_minmod;
		int flag_reconstruct = _REC_constant;
		int flag_gradient = _GRA_leastSquare;

	}
}

