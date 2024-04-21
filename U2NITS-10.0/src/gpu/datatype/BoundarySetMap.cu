#include "BoundarySetMap.h"

void GPU::BoundarySetMap::cuda_memcpy(BoundarySetMap* dist, const BoundarySetMap* src, cudaMemcpyKind kind) {
	integer num = dist->size;
	cudaMemcpy(dist->type, src->type, num * sizeof(integer), kind);
}
