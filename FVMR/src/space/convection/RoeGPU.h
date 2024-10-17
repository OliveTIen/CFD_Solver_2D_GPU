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
				// Harten-Yee型熵修正 目的是让特征值不要太接近0，导致非物理解
				// eig-特征值，lim-限制器，epsilon-防止分母为0
				// 当eig小于eig_lim时，增大eig，使其远离0
				// 基本不等式a^2+b^2>=2ab，因此(e*e+l*l)/(2*l)>=(2*e*l)/(2*l)= e

				// 原式： 
				if (eig < eig_lim) {
					eig = 0.5 * (eig * eig + eig_lim * eig_lim) / (eig_lim + epsilon);
				}

				//// 待验证？
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