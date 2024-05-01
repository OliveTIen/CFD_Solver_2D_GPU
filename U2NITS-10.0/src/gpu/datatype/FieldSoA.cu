#include "FieldSoA.h"
#include "../Env.h"

void GPU::ElementFieldSoA::cuda_memcpy(ElementFieldSoA* dist, const ElementFieldSoA* src, cudaMemcpyKind kind) {
	myint _num = dist->num;
	for (int i = 0; i < 4; i++) {
		cudaMemcpy(dist->U[i], src->U[i], _num * sizeof(myfloat), kind);
		cudaMemcpy(dist->Ux[i], src->Ux[i], _num * sizeof(myfloat), kind);
		cudaMemcpy(dist->Uy[i], src->Uy[i], _num * sizeof(myfloat), kind);
		cudaMemcpy(dist->Flux[i], src->Flux[i], _num * sizeof(myfloat), kind);
	}
}


void GPU::EdgeFieldSoA::cuda_memcpy(EdgeFieldSoA* dist, const EdgeFieldSoA* src, cudaMemcpyKind kind) {
	myint _num = dist->num_edge;
	for (int i = 0; i < 4; i++) {
		cudaMemcpy(dist->Flux[i], src->Flux[i], _num * sizeof(myfloat), kind);
	}
}
