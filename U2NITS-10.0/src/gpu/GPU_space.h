#ifndef GPU_SPACE_H
#define GPU_SPACE_H

#include "ElementDataPack.h"


namespace GPU {
	void calculateGradient(GPU::ElementDataPack& elementDataPack);

	__global__ void GPU::calculateGradientKernel(GPU::ElementDataPack& deviceDataPack);


	namespace Limiter {
		__device__ inline void Barth();
	}
}


#endif
