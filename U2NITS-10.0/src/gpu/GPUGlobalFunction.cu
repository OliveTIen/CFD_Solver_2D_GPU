#include "GPUGlobalFunction.h"
#include <iostream>
#include <ios>
#include "../output/LogWriter.h"
#include "../global/CExit.h"
#include "Env.h"

void GPU::catchCudaErrorAndExit() {
	cudaError_t cuda_error = cudaGetLastError();
	if (cuda_error != 0) {
		std::string e = "cudaError=" + std::to_string(cuda_error) + ", " + cudaGetErrorString(cuda_error);
		LogWriter::logAndPrint(e, LogWriter::Error, LogWriter::Error);
		CExit::saveAndExit(cuda_error);
	}



	//GPU::ErrorCode::ErrorCode errorCode = GPU::ErrorCode::noError;
	//cudaMemcpyFromSymbol(&errorCode, &GPU::globalErrorCode, sizeof(GPU::ErrorCode::ErrorCode), 0, cudaMemcpyDeviceToHost);
	//if (errorCode != GPU::ErrorCode::noError) {
	//	printf("ErrorCode: %d\n ErrorMessage: ", errorCode);
	//	printf(globalErrorMessage);
	//}
}

//void GPU::getGlobalErrorMessage(char* str, int size) {
//
//}
//
//__device__ GPU::ErrorCode::ErrorCode GPU::globalErrorCode = GPU::ErrorCode::noError;
//__device__ char GPU::globalErrorMessage[errorMessageLength];
