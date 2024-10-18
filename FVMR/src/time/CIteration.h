#ifndef CITERATION
#define CITERATION
#include "../gpu/datatype/DefineType.h"
#include "../solvers/GPUSolver2.h"

class CIteration {
public:
	static void iteration_host(myfloat& t, myfloat T, GPU::GPUSolver2* pSolver);///< iteration on CPU
	static void iteration_device(myfloat& t, myfloat T, GPU::GPUSolver2* pSolver);///< iteration on GPU

	// 更新host端的U_old和ruvp
	static void update_host_ruvp_Uold(const myfloat* const U[4], myfloat* ruvp[4], myfloat* U_old[4], myint num, myfloat gamma);
	static void update_device_ruvp_Uold(GPU::ElementFieldSoA& elementField_device, myfloat gamma);

private:
};

#endif // !CITERATION
