#include "BoundarySetMap.h"

void GPU::BoundarySetMap::cuda_memcpy(BoundarySetMap* dist, const BoundarySetMap* src, cudaMemcpyKind kind) {
	myint num = dist->size;
	cudaMemcpy(dist->type, src->type, num * sizeof(myint), kind);
}
