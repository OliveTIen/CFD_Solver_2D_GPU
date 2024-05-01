#ifndef _SPACE_RESTRICT_GPU_H_
#define _SPACE_RESTRICT_GPU_H_

#include "../../gpu/datatype/DataType.h"
#include "../../math/Common.h"

namespace GPU {
	namespace Space {
		namespace Restrict {

			inline __device__ bool outOfRange(myfloat* ruvp) {
				if (ruvp[0] < U2NITS::Math::Physics::RHO_MIN ||
					ruvp[0] > U2NITS::Math::Physics::RHO_MAX ||
					ruvp[3] < U2NITS::Math::Physics::P_MIN ||
					ruvp[3] > U2NITS::Math::Physics::P_MAX) {
					return true;
				}
				return false;
			}
		}
	}
}


#endif