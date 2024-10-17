#ifndef __EVOLVE_H__
#define __EVOLVE_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace U2NITS {
	namespace Time {

		// �����ƽ���ȫ��ʱ�䲽��
		void evolve_explicit_globaltimestep(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		// Runge-Kutta3��ȫ��ʱ�䲽��
		void evolve_rk3_globaltimestep(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		// ���㳣΢�ַ��̵��Ҷ���f=f(t,U).�����ع�
		void calculateFunctionF(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& ynp);

		// �����ƽ����ֲ�ʱ�䲽������Ҫdt����(λ��elementFieldVariable_dt_host.alphaC)
		void evolve_explicit_localtimestep(GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_host);

	}
}


#endif