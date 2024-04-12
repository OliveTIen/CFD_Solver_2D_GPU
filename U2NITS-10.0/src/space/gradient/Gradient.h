#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "../../gpu/datatype/DataType.h"

namespace U2NITS {
	namespace Space {
		namespace Gradient {

			void Gradient(GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::FieldSoA& elementField_host);
			void GradientLeastSquare_old(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host);
			void GradientLeastSquare_2(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::NodeSoA& node_host);
			void GradientGreenGauss_1(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host);
			void GradientGreenGauss_2(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host);

			//void LimiterBarth(GPU::ElementSoA& element_host, GPU::FieldSoA elementField_host);
		}
	}
}

#endif // !_GRADIENT_H_
