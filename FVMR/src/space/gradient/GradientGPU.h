#ifndef _CALCULATE_GRADIENT_2_CUH_
#define _CALCULATE_GRADIENT_2_CUH_

#include "../../gpu/datatype/DataType.h"

namespace GPU {
	void calculateGradient_old(GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device);
	namespace Space {
		namespace Gradient {
			void Gradient_2(GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device);
			void GradientLeastSquare(int block_size, int grid_size, GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device, GPU::NodeSoA& node_device);
			// �˺���������hostָ��(��������)�����޸�Ϊĳ��deviceָ��
			void GradientGreenGauss(int block_size, int grid_size, GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device, GPU::EdgeSoA& edge_device);
		}
	}

}

#endif