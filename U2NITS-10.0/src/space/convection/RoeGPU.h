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