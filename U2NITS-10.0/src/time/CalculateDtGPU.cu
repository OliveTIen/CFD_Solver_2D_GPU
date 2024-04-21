#include "CalculateDtGPU.h"
#include "../output/LogWriter.h"

real GPU::Time::calculateGlobalDt(real currentPhysicalTime, real maxPhysicalTime, real gamma, real Re, real Pr, real CFL, real Rcpcv, GPU::ElementSoA& element, GPU::EdgeSoA& edge, real* element_vruvp[4]) {
	/*
	涉及到规约操作
	*/
	LogWriter::logAndPrintError("unimplemented @GPU::Time::calculateGlobalDt\n");
	exit(-1);
}