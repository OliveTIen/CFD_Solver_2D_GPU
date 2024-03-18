#ifndef _CALCULATE_FLUX_2_CUH_
#define _CALCULATE_FLUX_2_CUH_

#include "../gpu/datatype/Datatype.h"
#include <map>
// __host__: 主机函数，由主机函数调用。是默认值
// __global__: 设备函数，由主机函数调用，且调用时要指明<<<>>>
// __device__: 设备函数，由设备函数调用

namespace GPU {
	namespace Space {
		// 未完成
		namespace Flux {
			void calculateFluxDevice(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DGlobalPara& infPara_device);

			void resetElementFlux(FieldSoA& elementField_device);
			void calculateFlux(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DGlobalPara& infPara_device);

			__global__ void resetElementFluxKernel(FieldSoA& elementField_device);
			__global__ void calculateFluxKernel(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DGlobalPara& infPara_device);
			__device__ void getEdgeFlux_wallNonViscous(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, long idx, REAL* flux, DGlobalPara& infPara_device);
			__device__ void getEdgeFlux_wall_adiabat(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, long idx, REAL* flux, DGlobalPara& infPara_device);
			__device__ void getEdgeFlux_farField(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, long idx, REAL* flux, DGlobalPara& infPara_device);
			__device__ void modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, DGlobalPara& infPara_device);
			__device__ void getEdgeFlux_inner_and_periodic(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, long idx, REAL* flux, DGlobalPara& infPara_device);
			__device__ void getUByXYandElementID(ElementSoA& element_device, FieldSoA& elementField_device, REAL x, REAL y, int elementID, REAL* U);

			__device__ void RiemannSolve(
				const REAL* UL, const REAL* UR, const REAL nx, const REAL ny,
				const REAL length, REAL* flux, DGlobalPara& para);
		}
	}
}


#endif