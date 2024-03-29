#include "GPUGlobalFunction.h"
#include <iostream>
#include <ios>
#include "../output/LogWriter.h"

void GPU::catchCudaErrorAndExit() {
	cudaError_t cuda_error = cudaGetLastError();
	if (cuda_error != 0) {
		std::string e = "cudaError=" + std::to_string(cuda_error) + ", " + cudaGetErrorString(cuda_error);
		LogWriter::logAndPrint(e, LogWriter::Error, LogWriter::Error);
		exit(cuda_error);
	}
}
