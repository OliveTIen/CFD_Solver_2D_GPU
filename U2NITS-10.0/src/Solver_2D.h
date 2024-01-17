// �����
// ����������������evolve(dt)�������ݻ�
// evolve(dt)�����ع����ݻ���ʽ���в�ͬ��ʽ������evolve_linear_explicit(dt)��evolve_linear_RK3(dt)
// ��Щ��������ͨ��(����calFlux())��Ȼ��ʱ���ƽ�
// calFlux()�ȵ���calEdgeFlux_Riemann(pEdge)���ÿ���ߵ�ͨ����Ȼ��Ӽ������൥Ԫ
// calEdgeFlux_Riemann(pEdge)���ݱ����ͣ�Ҳ�о������ʽ
//

#ifndef SOLVER_2D
#define SOLVER_2D

#include "./include/Eigen/Core" // �ú����ⱻU_2_F_lambdaʹ��
#include "./include/Eigen/Dense"
class Edge_2D;
class Element_T3;

class Solver_2D {
private:
	// static double RK3alpha[6];
	// static double RK5alpha[6];
	
public:
	// ��ת����flag=-1��ʾ�����
	// Eigen::Matrix4d get_matrixT(double nx, double ny, int flag = 1);
	void U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, double& lambda);

	// �ݻ�
	void evolve(double dt);
	// ��ʽ�ƽ�
	void evolve_explicit(double dt);
	// RK�ƽ�
	void evolve_RK3(double dt);


	//// [ʹ����]�����ֵͨ�� 
	void calFlux();//��ǰ
	void getEdgeFlux_inner(Edge_2D* pE, double* flux);
	void getEdgeFlux_wallNonViscous(Edge_2D* pE, double* flux);
	void getEdgeFlux_farfield(Edge_2D* pE, const double* ruvp_inf, double* flux);
	// ����Զ���߽����������ڳ�ֵ ����һ������ʹ��
	void modify_ruvpL_in_farfield(const double nx, const double ny, double* ruvp, const double* ruvp_inf);
	void getEdgeFlux_periodic(Edge_2D* pE, double* flux);
	/// �õ����Ӻ���
	

	//////ԭ�ȵ� ���������������ʦ˵�����д��ں�ɢ
	//void calFlux_old_LLF_1();
	//Eigen::Vector4d calEdgeFlux_LLF_old(Edge_2D* pEdge);
	//Eigen::Vector4d calEdgeFlux_LLF_wallNonViscous_old(Edge_2D* pEdge);
	//Eigen::Vector4d calEdgeFlux_LLF_farfield_old(Edge_2D* pEdge, const double* inf_ruvp);
	//Eigen::Vector4d calEdgeFlux_LLF_periodic_old(Edge_2D* pE);
	//Eigen::Vector4d calEdgeFlux_LLF_inner_old(Edge_2D* pE);
	//Eigen::Vector4d get_F_LLF_old(Eigen::Vector4d U_L, Eigen::Vector4d U_R, Edge_2D* pEdge);//[δʹ��]
	//Eigen::Vector4d get_Fv_old(Edge_2D* pEdge);//[δʹ��]ճ��ͨ��
	////��ʦ�Ĵ���
	//void calFlux_LLF_2();//��������ڱ߽� ��ʦ��fortran���� �������
	//void calFlux_Roe_2();//��������ڱ߽� ��ʦ��fortran���� ���ڷ�ɢ����
	void Compute_Deltaeig();//[δʹ��]Roe��ʽ������ֵ���㺯��

	// ������������ú������ڴ�������������ķ���ֵ����RiemannSolver����GlobalPara�ࡢ
	// LogWriter�����
	void RiemannSolve(const double* UL, const double* UR,
		const double nx, const double ny, const double length, double* flux,
		const int conservation_scheme);

	// [������]ճ��ͨ��
	void flux_viscous(Edge_2D* pEdge, double* flux_viscous);
};

#endif