#include "EdgeSoA.h"

void GPU::EdgeSoA::cuda_memcpy(EdgeSoA* dist, const EdgeSoA* src, cudaMemcpyKind kind) {
	myint num_edge = dist->num_edge;
	cudaMemcpy(dist->ID, src->ID, num_edge * sizeof(myint), kind);
	cudaMemcpy(dist->nodes[0], src->nodes[0], num_edge * sizeof(myint), kind);
	cudaMemcpy(dist->nodes[1], src->nodes[1], num_edge * sizeof(myint), kind);
	cudaMemcpy(dist->setID, src->setID, num_edge * sizeof(myint), kind);
	cudaMemcpy(dist->periodicPair, src->periodicPair, num_edge * sizeof(myint), kind);
	cudaMemcpy(dist->elementL, src->elementL, num_edge * sizeof(myint), kind);
	cudaMemcpy(dist->elementR, src->elementR, num_edge * sizeof(myint), kind);
	cudaMemcpy(dist->length, src->length, num_edge * sizeof(myfloat), kind);
	cudaMemcpy(dist->distanceOfElements, src->distanceOfElements, num_edge * sizeof(myfloat), kind);
	cudaMemcpy(dist->xy[0], src->xy[0], num_edge * sizeof(myfloat), kind);
	cudaMemcpy(dist->xy[1], src->xy[1], num_edge * sizeof(myfloat), kind);
	cudaMemcpy(dist->normal[0], src->normal[0], num_edge * sizeof(myfloat), kind);
	cudaMemcpy(dist->normal[1], src->normal[1], num_edge * sizeof(myfloat), kind);

}
