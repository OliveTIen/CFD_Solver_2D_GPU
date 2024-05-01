#ifndef PHYSICAL_KERNEL_H
#define PHYSICAL_KERNEL_H

#include "Common.h"

namespace U2NITS {
	namespace Math {

		// Ĭ�ϱ��Ϊglobal��global�ǲ��ǿ��Լȿ�����hostҲ������device���У�
		// ��U����ת��Ϊruvp����
		inline void U2ruvp_host(const myfloat U[4], myfloat ruvp[4], myfloat gamma) {
			// U:	rho,rho_u,rho_v,rho_E ע��E=e+0.5*V^2
			// ruvp: rho,u,v,p
			ruvp[0] = U[0];
			ruvp[1] = U[1] / U[0];
			ruvp[2] = U[2] / U[0];
			ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
		}

		inline void ruvp2U_host(const myfloat ruvp[4], myfloat U[4], myfloat gamma) {
			// ruvp: rho,u,v,p
			// U:	rho,rho_u,rho_v,rho_E
			U[0] = ruvp[0];
			U[1] = ruvp[0] * ruvp[1];
			U[2] = ruvp[0] * ruvp[2];
			U[3] = ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);//rhoE
		}

		// ��ά �غ���ת������
		inline void U2ruvwp_host_3d(const myfloat U[5], myfloat ruvwp[5], myfloat gamma) {

			/*
			// [��ɾ]ԭʼ�汾�����������Ƶ�����Ҫɾ�˸��Լ����鷳
			typedef myfloat myfloat;
			myfloat rho = U[0];
			myfloat rhou = U[1];
			myfloat rhov = U[2];
			myfloat rhow = U[3];
			myfloat rhoE = U[4];// ���������е��غ���
			myfloat u = rhou / rho;
			myfloat v = rhov / rho;
			myfloat w = rhow / rho;
			myfloat E = rhoE / rho;
			myfloat V2 = u * u + v * v + w * w;
			myfloat e = E - 0.5 * V2;
			// e = p/rho/(gamma-1)
			myfloat p = e * rho * (gamma - 1);
			*/

			// �򻯰汾
			ruvwp[0] = U[0];
			ruvwp[1] = U[1] / U[0];
			ruvwp[2] = U[2] / U[0];
			ruvwp[3] = U[3] / U[0];
			myfloat V2 = ruvwp[1] * ruvwp[1] + ruvwp[2] * ruvwp[2] + ruvwp[3] * ruvwp[3];
			ruvwp[4] = (U[4] / U[0] - 0.5 * V2) * U[0] * (gamma - 1);
		}

		// ĳ�����ϵ������
		inline myfloat getMach(const myfloat U[4], myfloat nx, myfloat ny, myfloat gamma) {
			myfloat ruvp[4]{};
			U2ruvp_host(U, ruvp, gamma);
			myfloat rho = ruvp[0];
			myfloat u = ruvp[1];
			myfloat v = ruvp[2];
			myfloat p = ruvp[3];
			myfloat un = u * nx + v * ny;
			myfloat a2 = gamma * p / rho;
			myfloat Ma2 = (u * u + v * v) / a2;
			return sqrt(Ma2);
		}

		//inline myfloat getPressureCoeffient() {
		//	// Cpp(I, J, K) = 2.0_8 * (PP(I, J, K) - PPF)   !!ѹ��ϵ��

		//}
	}
}


#endif