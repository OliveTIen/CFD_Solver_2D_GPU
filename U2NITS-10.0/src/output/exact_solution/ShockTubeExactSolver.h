#ifndef SHOCK_TUBE_EXACT_SOLVER_H
#define SHOCK_TUBE_EXACT_SOLVER_H
#include "../../gpu/datatype/DefineType.h"
class ShockTubeExactSolver {
public:
	void solveOneStep(myfloat* ruvp_inlet, myfloat* ruvp_outlet, myfloat* ruvp_output);
};

#endif // !SHOCK_TUBE_EXACT_SOLVER_H
