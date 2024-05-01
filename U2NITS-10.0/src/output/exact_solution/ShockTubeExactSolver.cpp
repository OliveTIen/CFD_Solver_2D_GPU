#include "ShockTubeExactSolver.h"
#include <iostream>
#include "../../global/GlobalPara.h"

void ShockTubeExactSolver::solveOneStep(myfloat* ruvp_inlet, myfloat* ruvp_outlet, myfloat* ruvp_output) {
	myfloat t{};
	myfloat x{};
	if (t == 0) {
		std::cout << "Error: t==0\n";
		return;
	}
	myfloat w = x / t;


}
