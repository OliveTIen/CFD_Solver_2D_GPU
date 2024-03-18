#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "../../gpu/datatype/DataType.h"

namespace U2NITS {
	namespace Space {
		namespace Gradient {
			void GradientLeastSquare(GPU::ElementSoA& element_host, GPU::FieldSoA elementField_host);
		}
	}
}

#endif // !_GRADIENT_H_
