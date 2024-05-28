#include "GPUGlobalFunction.h"
#include <iostream>
#include <sstream>
#include <ios>
#include "../output/LogWriter.h"
#include "../global/CExit.h"
#include "Env.h"

void GPU::catchCudaErrorAndExit() {
	cudaError_t cuda_error = cudaGetLastError();
	if (cuda_error != cudaSuccess) {
		std::string e = "cudaError=" + std::to_string(cuda_error) + ", " + cudaGetErrorString(cuda_error);
		LogWriter::logAndPrint(e, LogWriter::Error, LogWriter::Error);
		CExit::saveAndExit(cuda_error);
	}
}

void __getLastCudaError(const char* errorMessage, const char* file, const int line) {
	// 参见cuda samples项目中的helper_cuda.h
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err) {
		//// 原fprintf输出
		//fprintf(stderr,
		//	"%s(%i) : getLastCudaError() CUDA error :"
		//	" %s : (%d) %s.\n",
		//	file, line, errorMessage, static_cast<int>(err),
		//	cudaGetErrorString(err));
		// Log
		std::stringstream ss;
		ss << file << "(" << line << ")" << " : CUDA error : " << errorMessage << " : (" << static_cast<int>(err) << ") " << cudaGetErrorString(err) << ".\n";
		LogWriter::logAndPrintError(ss.str());
		exit(EXIT_FAILURE);
	}
}
