#ifndef SOLVER_DATA_GETTER_H
#define SOLVER_DATA_GETTER_H

namespace GPU {
	class GPUSolver2;
}

class SolverDataGetter {
public:
	// ��ȡSolverʵ����Ϊ������ϣ���ʡ���±���ʱ�䣬��SolverDataGetter��ӵ���
	// ����ָ֤��ǿգ�����֤��Դ������
	static GPU::GPUSolver2* getSolverInstance();
	//// ��ȡSolver����ָ��
	//static GPU::OutputNodeFieldSoA* getSolverOutputNodeFieldPointer();
};

#endif // !SOLVER_DATA_GETTER_H
