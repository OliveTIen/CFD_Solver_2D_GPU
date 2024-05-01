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
			void calculateFluxDevice(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, SDevicePara& para);


			void calculateFlux(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, SDevicePara& para);

			__global__ void calculateFluxKernel(ElementSoA element_device, ElementFieldSoA elementField_device, EdgeSoA edge_device, BoundarySetMap boundary_device, SDevicePara para);
			__device__ void getEdgeFlux_wallNonViscous_kernel(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, long idx, myfloat* flux, SDevicePara para);
			__device__ void getEdgeFlux_wall_adiabat_kernel(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, long idx, myfloat* flux, SDevicePara para);
			__device__ void getEdgeFlux_farField_kernel(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, long idx, myfloat* flux, SDevicePara para);
			__device__ void modify_ruvpL_farField_kernel(const myfloat nx, const myfloat ny, myfloat* ruvp, SDevicePara para);
			__device__ void getEdgeFlux_inner_and_periodic_kernel(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, long idx, myfloat* flux, SDevicePara para);
			__device__ void getUByXYandElementID_kernel(ElementSoA& element_device, ElementFieldSoA& elementField_device, myfloat x, myfloat y, int elementID, myfloat* U, SDevicePara para);

			__device__ void RiemannSolve_kernel(
				const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny,
				const myfloat length, myfloat* flux, SDevicePara para);

			void calculateFluxDevice_2(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, SDevicePara& para);
		}
	}
}


#endif