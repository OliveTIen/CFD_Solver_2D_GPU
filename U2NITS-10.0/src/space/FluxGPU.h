#ifndef _CALCULATE_FLUX_2_CUH_
#define _CALCULATE_FLUX_2_CUH_

#include "../gpu/datatype/Datatype.h"
#include <map>
// __host__: �����������������������á���Ĭ��ֵ
// __global__: �豸�������������������ã��ҵ���ʱҪָ��<<<>>>
// __device__: �豸���������豸��������

namespace GPU {
	namespace Space {
		// δ���
		namespace Flux {
			void calculateFluxDevice(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device);

			void resetElementFlux(FieldSoA& elementField_device);
			void calculateFlux(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device);
			__global__ void resetElementFluxKernel(FieldSoA& elementField_device);
			__global__ void calculateFluxKernel(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device);
			__device__ void getEdgeFlux_wallNonViscous(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux);
			__device__ void getEdgeFlux_wall_adiabat(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux);
			__device__ void getEdgeFlux_farField(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* ruvp_inf, REAL* flux);
			__device__ void modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, const REAL* ruvp_inf);
			__device__ void getEdgeFlux_inner_and_periodic(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* flux);
			__device__ void getUByXYandElementID(ElementSoA& element_host, FieldSoA& elementField_host, REAL x, REAL y, int elementID, REAL* U);

			__device__ void RiemannSolve(
				const REAL* UL, const REAL* UR, const REAL nx, const REAL ny,
				const REAL length, REAL* flux,
				const int scheme);
		}
	}
}


#endif