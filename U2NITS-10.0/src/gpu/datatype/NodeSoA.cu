#include "NodeSoA.h"

void GPU::NodeSoA::cuda_memcpy(NodeSoA* dist, const NodeSoA* src, cudaMemcpyKind kind) {
	int num_node = dist->num_node;
	cudaMemcpy(dist->ID, src->ID, num_node * sizeof(int), kind);
	cudaMemcpy(dist->xy[0], src->xy[0], num_node * sizeof(REAL), kind);
	cudaMemcpy(dist->xy[1], src->xy[1], num_node * sizeof(REAL), kind);
}
