#ifndef CALCULATE_DT_GPU_H
#define CALCULATE_DT_GPU_H
#include "../gpu/datatype/DataType.h"

namespace GPU {
	namespace Time {
		// 计算全局时间步长，存入elementFieldVariable_dt_device.dev_output[0]
		void get_global_dt_device(
			myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_device,
			GPU::ElementSoA& element_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device, int physicsModel_equation
		);
		// 计算局部时间步长，存入存入elementFieldVariable_dt_device.alphaC
		void get_local_dt_device(
			myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_device,
			GPU::ElementSoA& element_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device, int physicsModel_equation
		);
	}
}

#endif // !CALCULATE_DT_GPU_H
