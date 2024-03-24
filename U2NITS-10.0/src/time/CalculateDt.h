#ifndef _CALCULATE_DT_H_
#define _CALCULATE_DT_H_

#include "../space/FluxGPU.h"
#include "../math/Math.h"
#include "../solvers/GPUSolver2.h"
#include <iostream>

namespace GPU {
	// CPU�汾
	REAL calculateDt(
		REAL t, REAL T, REAL gamma, REAL Re, REAL prl, REAL CFL,
		REAL Rcpcv,// ���峣��R
		bool ifViscous, bool ifConstantViscous,
		GPU::GPUSolver2& gpuSolver);
}

#endif