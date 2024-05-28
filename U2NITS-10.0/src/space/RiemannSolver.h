#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H
#include "../gpu/datatype/DefineType.h"
/**
* 无粘黎曼求解器 先坐标变换，然后根据两侧U值计算界面无粘通量(准一维)
* 
* 为减少文件相互依赖，
* 1.conservation_scheme并没有直接用GlobalPara，
*   而是作为参数传递，这样可以缩短编译时间。
*   一旦GlobalPara变动，所有包含GlobalPara.h的文件都会被重新编译
* 2.对于错误处理，并没有直接调用LogWriter输出，而是输出返回值，让
*   调用它的函数处理错误。这样可以减少对LogWriter的依赖
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
