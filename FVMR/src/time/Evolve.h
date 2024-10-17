#ifndef __EVOLVE_H__
#define __EVOLVE_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace U2NITS {
	namespace Time {

		// 单步推进，全局时间步长
		void evolve_explicit_globaltimestep(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		// Runge-Kutta3，全局时间步长
		void evolve_rk3_globaltimestep(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		// 计算常微分方程的右端项f=f(t,U).包含重构
		void calculateFunctionF(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& ynp);

		// 单步推进，局部时间步长。需要dt数组(位于elementFieldVariable_dt_host.alphaC)
		void evolve_explicit_localtimestep(GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_host);

	}
}


#endif