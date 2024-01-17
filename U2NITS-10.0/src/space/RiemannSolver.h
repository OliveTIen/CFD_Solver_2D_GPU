#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H
/**
* ��ճ��������� ������任��Ȼ���������Uֵ���������ճͨ��(׼һά)
* 
* Ϊ�����ļ��໥������
* 1.conservation_scheme��û��ֱ����GlobalPara��
*   ������Ϊ�������ݣ������������̱���ʱ�䡣
*   һ��GlobalPara�䶯�����а���GlobalPara.h���ļ����ᱻ���±���
* 2.���ڴ�������û��ֱ�ӵ���LogWriter����������������ֵ����
*   �������ĺ�����������������Լ��ٶ�LogWriter������
*/
class RiemannSolver {
public:
	static enum ReturnStatus {
		normal,
		invalid_solver_type,
		compute_error
	};

public:
	static int solve(const double* UL, const double* UR,
		const double nx, const double ny, const double length, double* flux,
		const int conservation_scheme);
	static int LocalLaxFriedrichs(const double* UL, const double* UR, 
		const double nx, const double ny, const double length, double* flux);
};


#endif // !RIEMANN_SOLVER_H
