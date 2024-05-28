#ifndef GPU_GLOBAL_FUNCTION_H
#define GPU_GLOBAL_FUNCTION_H
#include "Env.h"
#include <stdio.h>

namespace GPU {
	constexpr int MY_BLOCK_SIZE = 512;// 512

	void catchCudaErrorAndExit();

	namespace ErrorCode {
		enum ErrorCode {
			noError,
			periodicPairNotFound
		};
	}


	inline void check_cudaError(cudaError_t status, const char* msg) {
		if (status != cudaSuccess) {
			const char* errorStr = cudaGetErrorString(status);
			printf("%s:\n%s\nError Code: %d\n\n", msg, errorStr, status);
			exit(status); // bail out immediately before the rest of the train derails...
		}
	}

	// Returns the maximum number of supported threads per block on the current device. Note: This is a hardware-enforced limit.
	inline int get_max_threads_per_block() {
		int dev_num;
		int max_threads;
		cudaError_t status;

		status = cudaGetDevice(&dev_num);// 当前活动设备
		check_cudaError(status, "Error querying device number.");

		status = cudaDeviceGetAttribute(&max_threads, cudaDevAttrMaxThreadsPerBlock, dev_num);
		check_cudaError(status, "Error querying max block threads.");

		return max_threads;
	}

}

	// This will output the proper error string when calling cudaGetLastError
#define getLastCudaError(msg) __getLastCudaError(msg, __FILE__, __LINE__)
	void __getLastCudaError(const char* errorMessage, const char* file, const int line);

#endif