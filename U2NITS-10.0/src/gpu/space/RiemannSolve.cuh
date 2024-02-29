#ifndef __RIEMANNSOLVE_CUH__
#define __RIEMANNSOLVE_CUH__

#include "../datatype/Define.h"

namespace GPU {
	namespace RiemannSolve {
		void RiemannSolve(
			const REAL* UL, const REAL* UR, const REAL nx, const REAL ny,
			const REAL length, REAL* flux,
			const int scheme);

		void LocalLaxFriedrichs(
			const REAL* UL, const REAL* UR, const REAL nx, const REAL ny,
			const REAL length, REAL* flux);

	}
}

#endif // __RIEMANNSOLVE_CUH__