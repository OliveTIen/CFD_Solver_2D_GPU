#ifndef _SPACE_RESTRICT_H_
#define _SPACE_RESTRICT_H_

#include "../../gpu/datatype/DataType.h"
#include "../../math/Common.h"

namespace U2NITS {
	namespace Space {
		namespace Restrict {
			// 修正所有单元场变量。对于超出范围的数据，取邻居的平均值
			void modifyElementFieldU2d(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField);
			// 修正iElemnet单元场变量。对于超出范围的数据，取邻居的平均值
			void modifyElementFieldUKernel2d(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, int iElement);

			inline bool outOfRange(myfloat* ruvp) {
				if (ruvp[0]<Math::Physics::RHO_MIN ||
					ruvp[0]>Math::Physics::RHO_MAX ||
					ruvp[3]<Math::Physics::P_MIN ||
					ruvp[3]>Math::Physics::P_MAX) {
					return true;
				}
				return false;
			}

			// 限制rho和p，使其非负且不太大
			inline void restrictRhoAndP(myfloat* ruvp) {
				if (ruvp[0] > Math::Physics::RHO_MAX)ruvp[0] = Math::Physics::RHO_MAX;
				if (ruvp[0] < Math::Physics::RHO_MIN)ruvp[0] = Math::Physics::RHO_MIN;
				if (ruvp[3] > Math::Physics::P_MAX)ruvp[3] = Math::Physics::P_MAX;
				if (ruvp[3] < Math::Physics::P_MIN)ruvp[3] = Math::Physics::P_MIN;
			}

			// 限制rho
			inline void restrictRho(myfloat* ruvpOrU) {
				if (ruvpOrU[0] > Math::Physics::RHO_MAX)ruvpOrU[0] = Math::Physics::RHO_MAX;
				if (ruvpOrU[0] < Math::Physics::RHO_MIN)ruvpOrU[0] = Math::Physics::RHO_MIN;
			}

			// 修正U，使其满足物理规律
			inline void restrictU2d(myfloat* U) {
				/*
				rhoE=rho(e+0.5*V2)，其中e=Cv*t>0，因此rhoE>0.5*rho*V2=0.5*1/rho*((rhoU)^2+(rhoV)^2)

				*/
				restrictRho(U);
				// 0.5*1/rho*((rhoU)^2+(rhoV)^2)
				myfloat rhoE_lowest = Math::EPSILON + 0.5 / U[0] * (U[1] * U[1] + U[2] * U[2]);
				U[3] = (U[3] > rhoE_lowest) ? U[3] : rhoE_lowest;
			}
			// 修正U
			inline void restrictU2d(myfloat& rho, myfloat& rhou, myfloat& rhov, myfloat& rhoE) {
				if (rho > Math::Physics::RHO_MAX)rho = Math::Physics::RHO_MAX;
				if (rho < Math::Physics::RHO_MIN)rho = Math::Physics::RHO_MIN;
				myfloat rhoE_lowest = Math::EPSILON + 0.5 / rho * (rhou * rhou + rhov * rhov);
				rhoE = (rhoE > rhoE_lowest) ? rhoE : rhoE_lowest;
			}
		}
	}
}

#endif // !_SPACE_RESTRICT_H_
