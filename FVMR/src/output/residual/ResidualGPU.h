#ifndef OUTPUT_RESIDUAL_GPU_H
#define OUTPUT_RESIDUAL_GPU_H
#include "../../gpu/datatype/FieldSoA.h"
#include "EnumNormType.h"
namespace GPU {
	namespace Output {
		// Î´Íê³É
		void get_residual_device(const GPU::ElementFieldSoA& elementField_device, myfloat* residual_device, U2NITS::Output::NormType NORM_TYPE);
	}
}
#endif // !OUTPUT_RESIDUAL_GPU_H
