#ifndef CITERATION
#define CITERATION
#include "../gpu/datatype/DefineType.h"
#include "../solvers/GPUSolver2.h"

class CIteration {
public:
	static void iteration_host_2(real& t, real T, GPU::GPUSolver2* pSolver);
	static void iteration_useGPU(real& t, real T, GPU::GPUSolver2* pSolver);


	static void update_host_ruvp_Uold(real* U[4], real* ruvp[4], real* U_old[4], integer num, real gamma);

private:
	static void m_EvolveHost_2_addRK3(real& physicalTime, real maxPhysicalTime, GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host, GPU::NodeSoA& node_host, real* element_vruvp[4]);
};

#endif // !CITERATION
