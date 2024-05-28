#ifndef OUTPUT_RESIDUAL_H
#define OUTPUT_RESIDUAL_H
#include "../../gpu/datatype/FieldSoA.h"
#include "EnumNormType.h"

namespace U2NITS {
	namespace Output {
		void get_residual_host(const GPU::ElementFieldSoA& elementField_host, myfloat* residual_host, NormType NORM_TYPE);
	}
}

#endif