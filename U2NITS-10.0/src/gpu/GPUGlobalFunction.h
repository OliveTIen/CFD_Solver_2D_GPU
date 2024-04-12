#ifndef GPU_GLOBAL_FUNCTION_H
#define GPU_GLOBAL_FUNCTION_H

#include "datatype/Define.h"

namespace GPU {
	void catchCudaErrorAndExit();

	namespace ErrorCode {
		enum ErrorCode {
			noError,
			periodicPairNotFound
		};
	}

	//constexpr int errorMessageLength = 200;
	//extern __device__ ErrorCode::ErrorCode globalErrorCode;
	//extern __device__ char globalErrorMessage[errorMessageLength];
	//void getGlobalErrorMessage(char* str, int size);
}

#endif