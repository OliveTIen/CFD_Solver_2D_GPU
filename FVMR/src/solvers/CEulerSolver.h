
#ifndef _CEULERSOLVER_H_
#define _CEULERSOLVER_H_
#include "CSolver.h"

namespace U2NITS {
	class CEulerSolver : CSolver {
	private:
		
	public:
		virtual void initialize();
		virtual void solve();
		virtual void updateResidual();
		virtual void finalize();
	};
}
#endif // !_CSOLVER_H_