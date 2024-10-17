#ifndef _CALCULATE_FLUX_2_CUH_
#define _CALCULATE_FLUX_2_CUH_

#include "../gpu/datatype/Datatype.h"
#include <map>
// __host__: 主机函数，由主机函数调用。是默认值
// __global__: 设备函数，由主机函数调用，且调用时要指明<<<>>>
// __device__: 设备函数，由设备函数调用

namespace GPU {
	namespace Space {
		// 未完成
		namespace Flux {
			void calculateFluxDevice_2(ElementSoA& element_device, EdgeSoA& edge_device, ElementFieldSoA& elementField_device);
		}
	}
}


#endif