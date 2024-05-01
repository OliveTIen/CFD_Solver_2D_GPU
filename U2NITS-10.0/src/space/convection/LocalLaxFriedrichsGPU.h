#ifndef LOCAL_LAX_FRIEDRICHS_GPU_H
#define LOCAL_LAX_FRIEDRICHS_GPU_H

#include "../../gpu/datatype/DefineType.h"
#include "../../gpu/Env.h"

namespace GPU {
	namespace Space {
		namespace Convection {
			__device__ void LocalLaxFriedrichs2d(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny,
				const myfloat length, myfloat* flux, myfloat gamma);
		}
	}
}

#endif // !LOCAL_LAX_FRIEDRICHS_H