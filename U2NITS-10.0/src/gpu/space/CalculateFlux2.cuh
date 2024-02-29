#ifndef _CALCULATE_FLUX_2_CUH_
#define _CALCULATE_FLUX_2_CUH_

#include "../dataType/ElementSoA.h"
#include "../dataType/EdgeSoA.h"
#include "../dataType/FieldSoA.h"
#include <map>

namespace GPU {
	// Î´Íê³É
	void calculateFlux2(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device);
	__global__ void calculateFluxKernel2(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device);

	void calculateFluxHost(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);

	REAL calculateLambdaFlux(REAL edgeU[4], REAL edgeN[2], REAL gamma, REAL length);

	namespace Flux {
		void getEdgeFlux_wallNonViscous(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux);
		void getEdgeFlux_wall_adiabat(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux);
		void getEdgeFlux_farField(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* ruvp_inf, REAL* flux);
		void modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, const REAL* ruvp_inf);
		void getEdgeFlux_inner_and_periodic(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair, long idx, REAL* flux);
		void getUByXYandElementID(ElementSoA& element_host, FieldSoA& elementField_host, REAL x, REAL y, int elementID, REAL* U);
	}
}

#endif