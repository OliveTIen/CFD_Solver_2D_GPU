#ifndef _SPACE_RESTRICT_GPU_H_
#define _SPACE_RESTRICT_GPU_H_

#include "../../gpu/datatype/DataType.h"
#include "../../math/CommonValue.h"

namespace GPU {
	namespace Space {
		namespace Restrict {

			inline __device__ bool outOfRange_device(myfloat* ruvp) {
				if (ruvp[0] < U2NITS::Math::Physics::RHO_MIN ||
					ruvp[0] > U2NITS::Math::Physics::RHO_MAX ||
					ruvp[3] < U2NITS::Math::Physics::P_MIN ||
					ruvp[3] > U2NITS::Math::Physics::P_MAX) {
					return true;
				}
				return false;
			}

			// 修正所有单元场变量。对于超出范围的数据，取邻居的平均值
			void modifyElementFieldU2d_device(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, const myfloat gamma);

		}
	}
}


#endif