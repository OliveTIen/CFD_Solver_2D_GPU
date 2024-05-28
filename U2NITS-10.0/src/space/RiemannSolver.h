#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H
#include "../gpu/datatype/DefineType.h"
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
	static int solve(const myfloat* UL, const myfloat* UR,
		const myfloat nx, const myfloat ny, const myfloat length, myfloat* flux,
		const int conservation_scheme);
	static int LocalLaxFriedrichs(const myfloat* UL, const myfloat* UR, 
		const myfloat nx, const myfloat ny, const myfloat length, myfloat* flux);
};


#endif // !RIEMANN_SOLVER_H
