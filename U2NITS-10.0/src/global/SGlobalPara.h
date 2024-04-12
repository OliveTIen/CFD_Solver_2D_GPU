#ifndef SGLOBALPARA_H
#define SGLOBALPARA_H
#include <string>
#include "Constexpr.h"
namespace SGlobalPara {


	struct SGlobalPara {
		/* struct type */

		struct Constant {
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

			void initialize(double _T0, double _p0, double _c0, double _gamma, double _epsilon, double _Re, double _Pr, double _mu) {
				T0 = _T0;
				p0 = _p0;
				c0 = _c0;
				gamma = _gamma;
				epsilon = _epsilon;
				Re = _Re;
				Pr = _Pr;
				mu = _mu;
			}
		};


		struct Basic {
			int dimension = 2;
			bool _continue = 1;
			std::string filename = "default_project";
			std::string meshFileType = "inp";
			bool useGPU = false;

		};


		struct Space {
			int flag_reconstruct = _REC_constant;
			int flag_gradient = _GRA_leastSquare;
		};

		struct Time {
			double CFL = 0.6;
			double T = 0.01;
			double residual = 1e-7;
			int time_advance = _EVO_explicit;

			// 以下数据从ContinueFile中读取，不是Toml
			double t_previous = 0;
			int istep_previous = 0;
		};

		struct PhysicsModel {
			int equation = 1;//"equation:1-Eluer,2-NS": 1
		};

		struct SingleBC_2D {
			bool use_ruvp = 0;
			double Ma = 0.8;
			double AOA = 1.25;//迎角
			double ruvp[4]{ 1,1,0,1 };
		};

		struct BoundaryCondition_2D {
			SingleBC_2D inlet;
			SingleBC_2D outlet;
			SingleBC_2D inf;
		};

		struct InitialCondition {
			int type = 1;
		};

		struct Output {
			int step_per_print = 50;
			int step_per_output_field = 50;
			int step_per_output_hist = 50;
			bool output_var_ruvp[4] = { 1,1,1,1 };
			int autosaveFileNum = 3;
			int maxIteration = 5000;
		};

		struct InviscidFluxMethod {
			int flux_conservation_scheme = _SOL_Roe;// 黎曼求解器
			int flux_limiter = _LIM_minmod;
		};

		/* members */
		Constant constant;
		Basic basic;
		Space space;
		Time time;
		PhysicsModel physicsModel;
		BoundaryCondition_2D boundaryCondition_2D;
		InitialCondition initialCondition;
		Output output;
		InviscidFluxMethod inviscidFluxMethod;
	};
}


#endif // !SGLOBALPARA_H
