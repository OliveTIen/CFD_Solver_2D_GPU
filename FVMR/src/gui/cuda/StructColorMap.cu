#include "StructColorMap.h"
#include "../../gpu/GPUGlobalFunction.h"

void GPU::ColorMap::cuda_memcpy(ColorMap* dist, const ColorMap* src, cudaMemcpyKind kind) {
	int num = src->num_control_point;
	cudaMemcpy(dist->data, src->data, num * sizeof(float4), kind);
}
