#include "FieldSoA.h"

void GPU::FieldSoA::cuda_memcpy(FieldSoA* dist, const FieldSoA* src, cudaMemcpyKind kind) {
	int _num = dist->num;
	for (int i = 0; i < 4; i++) {
		cudaMemcpy(dist->U[i], src->U[i], _num * sizeof(REAL), kind);
		cudaMemcpy(dist->Ux[i], src->Ux[i], _num * sizeof(REAL), kind);
		cudaMemcpy(dist->Uy[i], src->Uy[i], _num * sizeof(REAL), kind);
		cudaMemcpy(dist->Flux[i], src->Flux[i], _num * sizeof(REAL), kind);
	}
}

void GPU::OutputNodeFieldSoA::cuda_memcpy(OutputNodeFieldSoA* dist, const OutputNodeFieldSoA* src, cudaMemcpyKind kind) {
	int _num = dist->num_node;
	for (int i = 0; i < 4; i++) {
		cudaMemcpy(dist->ruvp[i], src->ruvp[i], _num * sizeof(REAL), kind);
	}
}

void GPU::EdgeFieldSoA::cuda_memcpy(EdgeFieldSoA* dist, const EdgeFieldSoA* src, cudaMemcpyKind kind) {
	int _num = dist->num_edge;
	for (int i = 0; i < 4; i++) {
		cudaMemcpy(dist->numeralFlux[i], src->numeralFlux[i], _num * sizeof(REAL), kind);
	}
}
