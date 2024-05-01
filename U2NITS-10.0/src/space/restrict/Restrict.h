#ifndef _SPACE_RESTRICT_H_
#define _SPACE_RESTRICT_H_

#include "../../gpu/datatype/DataType.h"
#include "../../math/Common.h"

namespace U2NITS {
	namespace Space {
		namespace Restrict {
			// �������е�Ԫ�����������ڳ�����Χ�����ݣ�ȡ�ھӵ�ƽ��ֵ
			void modifyElementFieldU2d(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField);
			// ����iElemnet��Ԫ�����������ڳ�����Χ�����ݣ�ȡ�ھӵ�ƽ��ֵ
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

			// ����rho��p��ʹ��Ǹ��Ҳ�̫��
			inline void restrictRhoAndP(myfloat* ruvp) {
				if (ruvp[0] > Math::Physics::RHO_MAX)ruvp[0] = Math::Physics::RHO_MAX;
				if (ruvp[0] < Math::Physics::RHO_MIN)ruvp[0] = Math::Physics::RHO_MIN;
				if (ruvp[3] > Math::Physics::P_MAX)ruvp[3] = Math::Physics::P_MAX;
				if (ruvp[3] < Math::Physics::P_MIN)ruvp[3] = Math::Physics::P_MIN;
			}

			// ����rho
			inline void restrictRho(myfloat* ruvpOrU) {
				if (ruvpOrU[0] > Math::Physics::RHO_MAX)ruvpOrU[0] = Math::Physics::RHO_MAX;
				if (ruvpOrU[0] < Math::Physics::RHO_MIN)ruvpOrU[0] = Math::Physics::RHO_MIN;
			}

			// ����U��ʹ�������������
			inline void restrictU2d(myfloat* U) {
				/*
				rhoE=rho(e+0.5*V2)������e=Cv*t>0�����rhoE>0.5*rho*V2=0.5*1/rho*((rhoU)^2+(rhoV)^2)

				*/
				restrictRho(U);
				// 0.5*1/rho*((rhoU)^2+(rhoV)^2)
				myfloat rhoE_lowest = Math::EPSILON + 0.5 / U[0] * (U[1] * U[1] + U[2] * U[2]);
				U[3] = (U[3] > rhoE_lowest) ? U[3] : rhoE_lowest;
			}
			// ����U
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
