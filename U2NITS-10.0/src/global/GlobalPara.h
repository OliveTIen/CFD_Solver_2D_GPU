#ifndef _CONSTEXPRS_H_
#define _CONSTEXPRS_H_

constexpr auto _NULL = 0;
//输出数据 1-999
constexpr auto _OUT_u = 1;
constexpr auto _OUT_v = 2;
constexpr auto _OUT_rho = 4;
constexpr auto _OUT_p = 8;
constexpr auto _OUT_T = 16;
constexpr auto _OUT_Unon_all = 101;
//方程类型 1001-1999
constexpr auto _EQ_convective = 1001;//对流扩散方程
constexpr auto _EQ_euler = 1002;//欧拉方程
constexpr auto _EQ_NS = 1003;
//求解器类型 2001-2999
constexpr auto _SOL_LocalLaxFriedrichs = 2001;//黎曼求解器 Local Lax-Friedrichs;
constexpr auto _SOL_Roe = 2002;
//重构方式 3001-3099
constexpr auto _REC_constant = 3001;//常量重构
constexpr auto _REC_linear = 3002;  //线性重构
//时间推进方式 3101-3199 evolution
constexpr auto _EVO_explicit = 3101;//显式时间推进
constexpr auto _EVO_rk3 = 3103;//RK3
constexpr auto _EVO_rk5 = 3105;
//限制器 4001-4999 limiter
constexpr auto _LIM_minmod = 4001;
//单元类型 5001-5999
constexpr auto _ELE_D2 = 5001;//2节点线段元
constexpr auto _ELE_T3 = 5101;//3节点三角元
constexpr auto _ELE_Q4 = 5102;//4节点四边元
//边界类型 6001-6999
constexpr auto _BC_wall_nonViscous = 6001;//固壁，无粘(有滑移)
constexpr auto _BC_wall_adiabat = 6002;//固壁，无滑移，绝热
constexpr auto _BC_wall_isothermal = 6003;//固壁，无滑移，等温
constexpr auto _BC_inlet = 6011;//入口
constexpr auto _BC_outlet = 6012;//出口
constexpr auto _BC_inf = 6013;//远场边界
constexpr auto _BC_symmetry = 6021;//对称边界
constexpr auto _BC_periodic_0 = 6100;//周期边界
constexpr auto _BC_periodic_1 = 6101;
constexpr auto _BC_periodic_2 = 6102;
constexpr auto _BC_periodic_3 = 6103;
constexpr auto _BC_periodic_4 = 6104;
constexpr auto _BC_periodic_5 = 6105;
constexpr auto _BC_periodic_6 = 6106;
constexpr auto _BC_periodic_7 = 6107;
constexpr auto _BC_periodic_8 = 6108;
constexpr auto _BC_periodic_9 = 6109;


#endif

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