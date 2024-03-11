/*!
 * \file roe.h
 * \brief Implementations of Roe-type schemes.
 * \author tgl
 *
 * 关于目录和文件的命名：
 * UNITs中，Roe位于Space/Convect_Roe.f90
 * FUN3D中，Roe位于libres/jacobian_gen.f90
 * SU2中，Roe位于src/numerics/flow/convection/roe.cpp
 * 
 */

#ifndef _ROE_H_
#define _ROE_H_
#include "../../gpu/datatype/Define.h"

namespace U2NITS {
	namespace Space {
		// 未使用
		void EigenValueAndVector4x4(REAL mat[4][4], REAL eigens[4], REAL R[4][4]);
		// 未使用
		void JacobiMethod(REAL mat[4][4], REAL eigens[4], REAL R[4][4]);
		// 未使用
		void RoeAverage(REAL U1[4], REAL U2[4], REAL gamma);
		// 未使用
		void RoeAverageFUN3D(
			REAL rho1, REAL rho2, REAL& rho,
			REAL u1, REAL u2, REAL& u,
			REAL v1, REAL v2, REAL& v,
			REAL h1, REAL h2, REAL& h,
			REAL e1, REAL e2, REAL& e,
			REAL frac1, REAL frac2, REAL& frac,
			REAL beta1, REAL beta2, REAL& beta,
			REAL& c, REAL& c2, REAL& q2,
			REAL p1, REAL p2,
			REAL n_species, REAL n_energy,
			REAL turb1, REAL turb2, REAL& turb,
			bool if_turb1);

		void ConvectRoeCommon3d(
			const REAL UL[5], const REAL UR[5], const REAL faceNormal[3],const REAL faceArea, REAL faceFlux[5],
			bool bDynamicMesh, REAL dynamicMeshValue, REAL gamma);
		// using by ConvectRoeCommon3d
		void RoeDissapationTerm3d(
			REAL gamma,
			const REAL UL[5], const REAL UR[5],
			REAL ruvwpL[5], REAL ruvwpR[5],
			const REAL faceNormal[3], REAL faceArea,
			bool bDynamicMesh, REAL dynamicMeshValue,
			bool bEntropyFix, REAL KEntropyFix[3], REAL kp,
			REAL drRoe[5]
		);

		void ConvectRoeCommon2d(
			const REAL UL[4], const REAL UR[4], const REAL faceNormal[2], const REAL faceArea, REAL faceFlux[4], REAL gamma
		);
		void RoeDissapationTerm2d(
			REAL gamma,
			REAL ruvpL[4], REAL ruvpR[4],
			const REAL faceNormal[2], REAL faceArea,
			bool bEntropyFix, REAL KEntropyFix[3], REAL p_sensor,
			REAL drRoe[4]
		);
		inline real PShockWaveSensor() { return 0.5; }


		// 未使用
		void GetRoeMatrix3d(
			REAL gamma, 
			REAL ruvwpL[5], REAL ruvwpR[5],
			REAL faceVector[3], REAL faceArea,
			bool bDynamicMesh, REAL dynamicMeshValue,
			bool bEntropyFix, REAL KEntropyFix[3],REAL kp,
			REAL roeMatrix[5][5]
		);

		inline void EigenEntropyFix_HartenYee(REAL& eig, REAL eig_lim, REAL epsilon) {
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

		inline REAL AdaptiveFunctionCoeffient() { return 1.0; }
	}
	namespace Matrix {
		// mat1: ixj, mat2: jxk, res: ixk
		inline void mul_ixj_jxk(int i, int j, int k, REAL* mat1, REAL* mat2, REAL* res) {
			for (int x = 0; x < i; x++) {
				for (int y = 0; y < k; y++) {
					*(res + x * k + y) = 0.0;
					for (int z = 0; z < j; z++) {
						// res[x][y] += mat1[x][z] * mat2[z][y];
						*(res + x * k + y) += *(mat1 + x * j + z) * *(mat2 + z * k + y);
					}
				}
			}
		}

	}
}


#endif