#include "NodeSoA.h"

void GPU::NodeSoA::cuda_memcpy(NodeSoA* dst, const NodeSoA* src, cudaMemcpyKind kind) {
	myint num_node = dst->num_node;
	cudaMemcpy(dst->ID, src->ID, num_node * sizeof(myint), kind);
	cudaMemcpy(dst->neighbor_element_count, src->neighbor_element_count, num_node * sizeof(myint), kind);
	cudaMemcpy(dst->xy[0], src->xy[0], num_node * sizeof(myfloat), kind);
	cudaMemcpy(dst->xy[1], src->xy[1], num_node * sizeof(myfloat), kind);
}
