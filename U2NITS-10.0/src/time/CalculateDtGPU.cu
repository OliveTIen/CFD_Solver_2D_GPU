#include "CalculateDtGPU.h"
#include "../output/LogWriter.h"

myfloat GPU::Time::calculateGlobalDt(myfloat currentPhysicalTime, myfloat maxPhysicalTime, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]) {
	/*
	涉及到规约操作
	*/
	LogWriter::logAndPrintError("unimplemented @GPU::Time::calculateGlobalDt\n");
	exit(-1);
}