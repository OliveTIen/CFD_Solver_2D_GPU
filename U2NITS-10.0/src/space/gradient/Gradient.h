#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "../../gpu/datatype/DataType.h"

namespace U2NITS {
	namespace Space {
		namespace Gradient {
			void Gradient(GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		}
	}
}

#endif // !_GRADIENT_H_
