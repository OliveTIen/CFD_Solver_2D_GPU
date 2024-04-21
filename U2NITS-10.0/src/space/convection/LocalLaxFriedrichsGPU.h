#ifndef LOCAL_LAX_FRIEDRICHS_GPU_H
#define LOCAL_LAX_FRIEDRICHS_GPU_H

#include "../../gpu/datatype/DefineType.h"

namespace GPU {
	namespace Space {
		namespace Convection {
			__device__ void LocalLaxFriedrichs2d(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny,
				const REAL length, REAL* flux, REAL gamma);
		}
	}
}

#endif // !LOCAL_LAX_FRIEDRICHS_H