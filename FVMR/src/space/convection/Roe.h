/*!
 * \file roe.h
 * \brief Implementations of Roe-type schemes.
 * \author tgl
 *
 * ����Ŀ¼���ļ���������
 * UNITs�У�Roeλ��Space/Convect_Roe.f90
 * FUN3D�У�Roeλ��libres/jacobian_gen.f90
 * SU2�У�Roeλ��src/numerics/flow/convection/roe.cpp
 * 
 */

#ifndef _ROE_H_
#define _ROE_H_
#include "../../gpu/datatype/DefineType.h"

namespace U2NITS {
	namespace Space {
		// δʹ��
		void EigenValueAndVector4x4(myfloat mat[4][4], myfloat eigens[4], myfloat R[4][4]);
		// δʹ��
		void JacobiMethod(myfloat mat[4][4], myfloat eigens[4], myfloat R[4][4]);
		// δʹ��
		void RoeAverage(myfloat U1[4], myfloat U2[4], myfloat gamma);
		// δʹ��
		void RoeAverageFUN3D(
			myfloat rho1, myfloat rho2, myfloat& rho,
			myfloat u1, myfloat u2, myfloat& u,
			myfloat v1, myfloat v2, myfloat& v,
			myfloat h1, myfloat h2, myfloat& h,
			myfloat e1, myfloat e2, myfloat& e,
			myfloat frac1, myfloat frac2, myfloat& frac,
			myfloat beta1, myfloat beta2, myfloat& beta,
			myfloat& c, myfloat& c2, myfloat& q2,
			myfloat p1, myfloat p2,
			myfloat n_species, myfloat n_energy,
			myfloat turb1, myfloat turb2, myfloat& turb,
			bool if_turb1);

		void ConvectRoeCommon3d(
			const myfloat UL[5], const myfloat UR[5], const myfloat faceNormal[3], const myfloat faceArea, myfloat faceFlux[5],
			bool bDynamicMesh, myfloat dynamicMeshValue, myfloat gamma, myfloat rcpcv);
		// using by ConvectRoeCommon3d
		void RoeDissapationTerm3d(
			myfloat gamma,
			const myfloat UL[5], const myfloat UR[5],
			myfloat ruvwpL[5], myfloat ruvwpR[5],
			const myfloat faceNormal[3], myfloat faceArea,
			bool bDynamicMesh, myfloat dynamicMeshValue,
			bool bEntropyFix, myfloat KEntropyFix[3], myfloat kp,
			myfloat drRoe[5]
		);

		void ConvectRoeCommon2d(
			const myfloat UL[4], const myfloat UR[4], const myfloat faceNormal[2], const myfloat faceArea, myfloat faceFlux[4], myfloat gamma, myfloat rcpcv
		);
		void RoeDissapationTerm2d(
			myfloat gamma,
			myfloat ruvpL[4], myfloat ruvpR[4],
			const myfloat faceNormal[2], myfloat faceArea,
			bool bEntropyFix, myfloat KEntropyFix[3], myfloat p_sensor,
			myfloat drRoe[4]
		);
		inline myfloat PShockWaveSensor() { return 0.5; }


		// δʹ��
		void GetRoeMatrix3d(
			myfloat gamma,
			myfloat ruvwpL[5], myfloat ruvwpR[5],
			myfloat faceVector[3], myfloat faceArea,
			bool bDynamicMesh, myfloat dynamicMeshValue,
			bool bEntropyFix, myfloat KEntropyFix[3], myfloat kp,
			myfloat roeMatrix[5][5]
		);

		inline void EigenEntropyFix_HartenYee(myfloat& eig, myfloat eig_lim, myfloat epsilon) {
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

		inline myfloat AdaptiveFunctionCoeffient() { return 1.0; }
	}
	namespace Matrix {
		// mat1: ixj, mat2: jxk, res: ixk
		inline void mul_ixj_jxk(int i, int j, int k, myfloat* mat1, myfloat* mat2, myfloat* res) {
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