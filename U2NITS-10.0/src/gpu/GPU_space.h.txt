#ifndef GPU_SPACE_H
#define GPU_SPACE_H

#include "dataType/ElementDataPack.h"


namespace GPU {
	void calculateGradient(GPU::ElementDataPack& element_device);

	__global__ void GPU::calculateGradientKernel(GPU::ElementDataPack& element_device);

	void updateNeighborGradient(GPU::ElementDataPack& element_host, GPU::ElementAdjacent& adjacent);

	namespace Limiter {
		__device__ inline void Barth();
	}

	void calculateFlux(GPU::ElementDataPack& element_device);
}


#endif
