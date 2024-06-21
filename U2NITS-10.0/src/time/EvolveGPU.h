#ifndef __EVOLVE_GPU_H__
#define __EVOLVE_GPU_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace GPU {
	namespace Time {
		// 单步推进，全局时间步长
		void evolve_explicit_globaltimestep_device(myfloat dt, GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device);
		// 单步推进，局部时间步长。需要dt数组(位于elementFieldVariable_dt_device.alphaC)
		void evolve_explicit_localtimestep_device(GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_device);

	}
}



#endif