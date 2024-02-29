#ifndef _CALCULATE_DT_H_
#define _CALCULATE_DT_H_

#include "../datatype/ElementSoA.h"
#include "../space/CalculateFlux2.cuh"
#include "../math/PhysicalConvertKernel.h"
#include <iostream>
#include "../GPUSolver2.h"

namespace GPU {
	// CPU版本
	REAL calculateDt(
		REAL t, REAL T, REAL gamma, REAL Re, REAL prl, REAL CFL,
		REAL Rcpcv,// 气体常数R
		bool ifViscous, bool ifConstantViscous,
		GPU::GPUSolver2& gpuSolver);
}

#endif