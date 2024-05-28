#ifndef CALCULATE_DT_GPU_H
#define CALCULATE_DT_GPU_H
#include "../gpu/datatype/DataType.h"

namespace GPU {
	namespace Time {
		// 计算全局时间步长，并存入elementFieldVariable_dt_device.dev_output[0]
		void get_global_dt_device(
			myfloat t, myfloat T, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_device,
			GPU::ElementSoA& element_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device, int physicsModel_equation
		);
	}
}

#endif // !CALCULATE_DT_GPU_H
