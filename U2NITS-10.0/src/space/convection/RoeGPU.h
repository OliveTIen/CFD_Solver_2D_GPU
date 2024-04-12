#ifndef _ROE_GPU_H_
#define _ROE_GPU_H_
#include "../../gpu/datatype/Datatype.h"


namespace GPU {
	namespace Space {
		namespace Convection {
			__device__ void ConvectRoeCommon2d(
				const REAL UL[4], const REAL UR[4], const REAL faceNormal[2], const REAL faceArea, REAL faceFlux[4], SDevicePara para
			);
			__device__ void RoeDissapationTerm2d(
				REAL gamma,
				REAL ruvpL[4], REAL ruvpR[4],
				const REAL faceNormal[2], REAL faceArea,
				bool bEntropyFix, REAL KEntropyFix[3], REAL p_sensor,
				REAL drRoe[4]
			);
			__device__ inline real PShockWaveSensor() { return 0.5; }

			__device__ inline void EigenEntropyFix_HartenYee(REAL& eig, REAL eig_lim, REAL epsilon) {
				// Harten-Yee�������� Ŀ����������ֵ��Ҫ̫�ӽ�0�����·������
				// eig-����ֵ��lim-��������epsilon-��ֹ��ĸΪ0
				// ��eigС��eig_limʱ������eig��ʹ��Զ��0
				// ��������ʽa^2+b^2>=2ab�����(e*e+l*l)/(2*l)>=(2*e*l)/(2*l)= e

				// ԭʽ�� 
				if (eig < eig_lim) {
					eig = 0.5 * (eig * eig + eig_lim * eig_lim) / (eig_lim + epsilon);
				}

				//// ����֤��
				//if (eig_lim < 0)eig_lim = -eig_lim;
				//if (eig_lim == 0)eig_lim = epsilon;
				//if (eig < eig_lim) {
				//	eig = 0.5 * (eig * eig + eig_lim * eig_lim) / (eig_lim + epsilon);
				//}
			}
		}
	}
}

#endif