#ifndef SHOCK_TUBE_EXACT_SOLVER_H
#define SHOCK_TUBE_EXACT_SOLVER_H
#include "../../gpu/datatype/DefineType.h"
class ShockTubeExactSolver {
public:
	void solveOneStep(real* ruvp_inlet, real* ruvp_outlet, real* ruvp_output);
};

#endif // !SHOCK_TUBE_EXACT_SOLVER_H
