#include "ElementSoA.h"



void GPU::ElementSoA::cuda_memcpy(ElementSoA* dist, const ElementSoA* src, cudaMemcpyKind kind) {
	integer num_element = dist->num_element;
	cudaMemcpy(dist->ID, src->ID, num_element * sizeof(integer), kind);
	cudaMemcpy(dist->xy[0], src->xy[0], num_element * sizeof(REAL), kind);
	cudaMemcpy(dist->xy[1], src->xy[1], num_element * sizeof(REAL), kind);
	cudaMemcpy(dist->volume, src->volume, num_element * sizeof(REAL), kind);
	for (integer i = 0; i < 4; i++) {
		cudaMemcpy(dist->nodes[i], src->nodes[i], num_element * sizeof(integer), kind);
		cudaMemcpy(dist->edges[i], src->edges[i], num_element * sizeof(integer), kind);
		cudaMemcpy(dist->neighbors[i], src->neighbors[i], num_element * sizeof(integer), kind);
		//cudaMemcpy(dist->U[i], src->U[i], num_element * sizeof(REAL), kind);
		//cudaMemcpy(dist->Ux[i], src->Ux[i], num_element * sizeof(REAL), kind);
		//cudaMemcpy(dist->Uy[i], src->Uy[i], num_element * sizeof(REAL), kind);
		//cudaMemcpy(dist->Flux[i], src->Flux[i], num_element * sizeof(REAL), kind);
	}
	
}
/*

	需要复制的变量：
		int num_element;

		int* ID;
		int* nodes[4];
		int* edges[4];
		int* neighbors[4];
		REAL* x;
		REAL* y;
		REAL* U[4];
		REAL* Ux[4];
		REAL* Uy[4];
		REAL* Flux[4];

	函数参照：
	int num_edge = src->num_edge;
	cudaMemcpy(dist->ID, src->ID, num_edge * sizeof(int), kind);
	cudaMemcpy(dist->nodes[0], src->nodes[0], num_edge * sizeof(int), kind);
	cudaMemcpy(dist->nodes[1], src->nodes[1], num_edge * sizeof(int), kind);
	cudaMemcpy(dist->setID, src->setID, num_edge * sizeof(int), kind);
	cudaMemcpy(dist->elementL, src->elementL, num_edge * sizeof(int), kind);
	cudaMemcpy(dist->elementR, src->elementR, num_edge * sizeof(int), kind);
	cudaMemcpy(dist->length, src->length, num_edge * sizeof(REAL), kind);
	cudaMemcpy(dist->distanceOfElements, src->distanceOfElements, num_edge * sizeof(REAL), kind);

*/

