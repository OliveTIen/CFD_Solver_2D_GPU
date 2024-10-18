#ifndef CALCULATE_DT_GPU_H
#define CALCULATE_DT_GPU_H
#include "../gpu/datatype/DataType.h"

namespace GPU {
	/** @brief namespace time	*/
	namespace Time {
		
		/** @brief get global time step. 
		* 计算全局时间步长，存入elementFieldVariable_dt_device.dev_output[0]
		* @details 先计算局部时间步长，然后取最小值<br>
		* 注意reduce_device中n必须为2的幂，用num_reduce而不是num_element
		* @param[in] gamma gamma
		* @param[in] Re Reynold number
		* @param[in] Pr Prantl number
		* @param[in] CFL Courant number
		* @param[in] Rcpcv specific heat ratio
		*
		*/
		void get_global_dt_device(
			myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_device,
			GPU::ElementSoA& element_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device, int physicsModel_equation
		);

		/** @brief get local time step.
		* 计算局部时间步长，存入elementFieldVariable_dt_device.alphaC
		* @param[in] gamma gamma
		* @param[in] Re Reynold number
		* @param[in] Pr Prantl number
		* @param[in] CFL Courant number
		* @param[in] Rcpcv specific heat ratio
		*
		*/
		void get_local_dt_device(
			myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_device,
			GPU::ElementSoA& element_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device, int physicsModel_equation
		);
	}
}

#endif // !CALCULATE_DT_GPU_H
