#ifndef SOLVER_DATA_GETTER_H
#define SOLVER_DATA_GETTER_H

namespace GPU {
	class GPUSolver2;
}

class SolverDataGetter {
public:
	// 获取Solver实例。为减少耦合，节省重新编译时间，用SolverDataGetter间接调用
	// 仅保证指针非空，不保证资源已申请
	static GPU::GPUSolver2* getSolverInstance();
	//// 获取Solver数据指针
	//static GPU::OutputNodeFieldSoA* getSolverOutputNodeFieldPointer();
};

#endif // !SOLVER_DATA_GETTER_H
