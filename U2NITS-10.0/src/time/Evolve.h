#ifndef __EVOLVE_H__
#define __EVOLVE_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace U2NITS {
	namespace Time {
		// 旧函数，已弃用
		void EvolveHost_1(myfloat dt, int flag, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);
		void EvolveExplicitHost(myfloat dt, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);

		// 非定常
		void evolve_unsteady_explicit(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		void evolve_unsteady_rk3(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		// 计算常微分方程的右端项f=f(t,U).包含重构
		void calculateFunctionF(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& ynp);
		void U_to_ruvp_proMax(const myfloat* U[4], myfloat* ruvp[4], int length, myfloat gamma);

		// 定常
		void evolve_steady_explicit_localTimeStep(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, myfloat* element_vruvp[4]);


	}
}


#endif