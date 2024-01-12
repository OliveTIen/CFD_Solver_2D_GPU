// 求解器
// 主程序调用求解器的evolve(dt)，进行演化
// evolve(dt)根据重构和演化方式，有不同形式，例如evolve_linear_explicit(dt)、evolve_linear_RK3(dt)
// 这些函数先求通量(调用calFlux())，然后时间推进
// calFlux()先调用calEdgeFlux_Riemann(pEdge)求解每个边的通量，然后加减到两侧单元
// calEdgeFlux_Riemann(pEdge)根据边类型，也有具体的形式
//

#ifndef SOLVER_2D
#define SOLVER_2D

#include "./include/Eigen/Core" // 该函数库被U_2_F_lambda使用
#include "./include/Eigen/Dense"
class Edge_2D;
class Element_T3;

class Solver_2D {
public:
	//static double RK3alpha[6];
	//static double RK5alpha[6];
	
public:
	//旋转矩阵。flag=-1表示逆矩阵
	//Eigen::Matrix4d get_matrixT(double nx, double ny, int flag = 1);
	void U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, double& lambda);

	//演化
	void evolve(double dt);
	//显式推进
	void evolve_explicit(double dt);
	//RK推进
	void evolve_RK3(double dt);


	////求解数值通量 new 将LocalLaxFriedrichs和Roe放进其中
	void calFlux();

	////[使用中]将LLF和Roe集成
	void calFlux_current();//当前
	void getEdgeFlux_inner(Edge_2D* pE, double* flux);
	void getEdgeFlux_wallNonViscous(Edge_2D* pE, double* flux);
	void getEdgeFlux_farfield(Edge_2D* pE, const double* ruvp_inf, double* flux);
	void getEdgeFlux_periodic(Edge_2D* pE, double* flux);
	///用到的子函数
	
	void LLF_new_current(const double* UL, const double* UR, const double nx, const double ny, const double length, double* flux);//功能：根据ULUR等参数，计算flux[0-3]
	void LLF_test_2(const double* UL, const double* UR, const double nx, const double ny, const double length, double* flux);//测试
	void LLF_test(const double* UL, const double* UR, const double nx, const double ny, const double length, double* flux);//测试
	void cal_ruvp_farfield_new(const double nx, const double ny, double* ruvp, const double* ruvp_inf);

	//////原先的 能算出来，但是老师说等熵涡存在耗散
	//void calFlux_old_LLF_1();
	//Eigen::Vector4d calEdgeFlux_LLF_old(Edge_2D* pEdge);
	//Eigen::Vector4d calEdgeFlux_LLF_wallNonViscous_old(Edge_2D* pEdge);
	//Eigen::Vector4d calEdgeFlux_LLF_farfield_old(Edge_2D* pEdge, const double* inf_ruvp);
	//Eigen::Vector4d calEdgeFlux_LLF_periodic_old(Edge_2D* pE);
	//Eigen::Vector4d calEdgeFlux_LLF_inner_old(Edge_2D* pE);
	//Eigen::Vector4d get_F_LLF_old(Eigen::Vector4d U_L, Eigen::Vector4d U_R, Edge_2D* pEdge);//[未使用]
	//Eigen::Vector4d get_Fv_old(Edge_2D* pEdge);//[未使用]粘性通量


	////老师的代码
	//void calFlux_LLF_2();//仅完成周期边界 老师的fortran代码 能算出来
	//void calFlux_Roe_2();//仅完成周期边界 老师的fortran代码 存在发散问题
	void Compute_Deltaeig();
};

#endif