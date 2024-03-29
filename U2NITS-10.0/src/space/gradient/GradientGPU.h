#ifndef _CALCULATE_GRADIENT_2_CUH_
#define _CALCULATE_GRADIENT_2_CUH_

#include "../../gpu/datatype/DataType.h"

namespace GPU {
	void calculateGradient_old(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device);

	__global__ void calculateGradientKernel_old(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device);

	namespace Space {
		namespace Gradient {
			void Gradient_2(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device, GPU::EdgeSoA& edge_device);
			void GradientLeastSquare(int block_size, int grid_size, GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device, GPU::NodeSoA& node_device);
			// 核函数不允许传host指针(包括引用)，需修改为某个device指针
			__global__ void GradientLeastSquareKernel(GPU::ElementSoA element_device, GPU::FieldSoA elementField_device, GPU::NodeSoA node_device);
			void GradientGreenGauss(int block_size, int grid_size, GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device, GPU::EdgeSoA& edge_device);
			__global__ void GradientGreenGaussKernel(GPU::ElementSoA element_device, GPU::FieldSoA elementField_device, GPU::EdgeSoA edge_device);
		}
	}

}

#endif