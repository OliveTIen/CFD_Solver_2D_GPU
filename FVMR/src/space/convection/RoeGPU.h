#ifndef _ROE_GPU_H_
#define _ROE_GPU_H_
#include "../../gpu/datatype/Datatype.h"


namespace GPU {
	namespace Space {
		namespace Convection {
			__device__ void ConvectRoeCommon2d(
				const myfloat UL[4], const myfloat UR[4], const myfloat faceNormal[2], const myfloat faceArea, myfloat faceFlux[4], myfloat gamma, myfloat rcpcv
			);
			__device__ void RoeDissapationTerm2d(
				myfloat gamma,
				myfloat ruvpL[4], myfloat ruvpR[4],
				const myfloat faceNormal[2], myfloat faceArea,
				bool bEntropyFix, myfloat KEntropyFix[3], myfloat p_sensor,
				myfloat drRoe[4]
			);
			__device__ inline myfloat PShockWaveSensor() { return 0.5; }

			__device__ inline void EigenEntropyFix_HartenYee(myfloat& eig, myfloat eig_lim, myfloat epsilon) {
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