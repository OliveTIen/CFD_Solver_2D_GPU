#include "ShockTubeExactSolver.h"
#include <iostream>
#include "../../global/GlobalPara.h"

void ShockTubeExactSolver::solveOneStep(real* ruvp_inlet, real* ruvp_outlet, real* ruvp_output) {
	real t{};
	real x{};
	if (t == 0) {
		std::cout << "Error: t==0\n";
		return;
	}
	real w = x / t;


}
