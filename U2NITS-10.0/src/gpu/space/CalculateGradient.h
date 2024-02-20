
#ifndef _CALCULATE_GRADIENT_H_
#define _CALCULATE_GRADIENT_H_

#include "../dataType/ElementDataPack.h"
#include "../dataType/ElementAdjacent.h"

// 报错：
// qualified name is not allowed in namespace member declaration
// 原来是函数声明前面误加了GPU::，

namespace GPU {
	void calculateGradient(GPU::ElementDataPack& element_device);

	__global__ void calculateGradientKernel(GPU::ElementDataPack& element_device);

	void updateNeighborGradient(GPU::ElementDataPack& element_host, GPU::ElementAdjacent& adjacent);

	namespace Limiter {
		__device__ inline void Barth();
	}



}

#endif