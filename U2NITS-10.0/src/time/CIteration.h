#ifndef CITERATION
#define CITERATION
#include "../gpu/datatype/DefineType.h"
#include "../solvers/GPUSolver2.h"

class CIteration {
public:
	static void iteration_host_2(myfloat& t, myfloat T, GPU::GPUSolver2* pSolver);
	static void iteration_useGPU(myfloat& t, myfloat T, GPU::GPUSolver2* pSolver);


	static void update_host_ruvp_Uold(myfloat* U[4], myfloat* ruvp[4], myfloat* U_old[4], myint num, myfloat gamma);

private:
	static void m_EvolveHost_2_addRK3(myfloat& physicalTime, myfloat maxPhysicalTime, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, GPU::NodeSoA& node_host, myfloat* element_vruvp[4]);
};

#endif // !CITERATION
