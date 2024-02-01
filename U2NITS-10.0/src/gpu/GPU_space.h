#ifndef GPU_SPACE_H
#define GPU_SPACE_H

#include "ElementDataPack.h"


namespace GPU {
	void calculateGradient(int block_size, int grid_size, GPU::ElementDataPack& elementDataPack);

	__global__ void GPU::calculateGradientKernel(GPU::ElementDataPack& deviceDataPack);
}


#endif
