#ifndef __EVOLVE_H__
#define __EVOLVE_H__

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace U2NITS {
	namespace Time {
		// �ɺ�����������
		void EvolveHost_1(REAL dt, int flag, GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);
		void EvolveExplicitHost(REAL dt, GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair);

		// �Ƕ���
		void evolveSingleStep(real dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::FieldSoA& elementField_host);
		void evolveRungeKutta3(real dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::FieldSoA& elementField_host);
		void calculateFunctionF(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::FieldSoA& ynp);// �����Ҷ�������ع�
		void U_to_ruvp_proMax(const real* U[4], real* ruvp[4], int length, real gamma);

		// ����
		void evolveSteadyLocalTimeStep(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::FieldSoA& elementField, real* element_vruvp[4]);


	}
}


#endif