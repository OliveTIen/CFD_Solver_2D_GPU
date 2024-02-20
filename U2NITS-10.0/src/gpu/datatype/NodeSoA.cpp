#include "NodeSoA.h"

void GPU::NodeSoA::cuda_memcpy(NodeSoA* dist, const NodeSoA* src, cudaMemcpyKind kind) {
	int num_node = dist->num_node;
	cudaMemcpy(dist->ID, src->ID, num_node * sizeof(int), kind);
	cudaMemcpy(dist->x, src->x, num_node * sizeof(REAL), kind);
	cudaMemcpy(dist->x, src->x, num_node * sizeof(REAL), kind);
}
