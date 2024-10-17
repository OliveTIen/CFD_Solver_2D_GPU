#include "Roe.h"
#include "../../math/Math.h"
#include "../../global/GlobalPara.h"
#include "../../output/LogWriter.h"
#include <cmath>
#include "../../global/CExit.h"

void U2NITS::Space::EigenValueAndVector4x4(myfloat mat[4][4], myfloat eigens[4], myfloat R[4][4]) {
	// ����4x4���������ֵ����������
	// ����ֱ�Ӳ���UNITs�ķ�����
}

void U2NITS::Space::JacobiMethod(myfloat mat[4][4], myfloat eigens[4], myfloat R[4][4]) {
	// ��R��ʼ��Ϊ��λ����
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				R[i][j] = 1.0;
			}
			else {
				R[i][j] = 0.0;
			}
		}
	}

	const int n = 4;
	const myfloat eps = 1.0e-10; // ��������
	myfloat max_off_diag; // ���ǶԽ�Ԫ��
	myfloat theta, t, c, s; // ��ת�ǶȺ���ת����
	int p, q;
	do {
		max_off_diag = 0.0;
		// Ѱ�����ǶԽ�Ԫ��
		for (int i = 0; i < n - 1; i++) {
			for (int j = i + 1; j < n; j++) {
				if (std::abs(mat[i][j]) > max_off_diag) {
					max_off_diag = std::abs(mat[i][j]);
					p = i;
					q = j;
				}
			}
		}

		if (max_off_diag < eps) {
			break; // �Ѿ��㹻�ӽ��Խǻ�
		}

		// ������ת�Ƕ�
		theta = 0.5 * std::atan2(2 * mat[p][q], mat[p][p] - mat[q][q]);// atan2���ط�λ��
		t = mat[p][q];
		c = std::cos(theta);
		s = std::sin(theta);

		// �޸ľ��������������
		mat[p][q] = 0.0;
		for (int i = 0; i < n; i++) {
			if (i != p && i != q) {
				myfloat temp = mat[i][p];
				mat[i][p] = temp * c - mat[i][q] * s;
				mat[p][i] = mat[i][p];
				mat[i][q] = mat[i][q] * c + temp * s;
				mat[q][i] = mat[i][q];
			}
			// ��������������
			myfloat new_R_ip = R[i][p] * c - R[i][q] * s;
			myfloat new_R_iq = R[i][q] * c + R[i][p] * s;
			R[i][p] = new_R_ip;
			R[i][q] = new_R_iq;
		}
		mat[p][p] = mat[p][p] * c * c - 2 * t * c * s + mat[q][q] * s * s;
		mat[q][q] = mat[q][q] * c * c + 2 * t * c * s + mat[p][p] * s * s;
		mat[p][p] = mat[q][q]; // ��������СԪ��
		mat[q][q] = c * c * mat[q][q] + 2 * t * c * s + s * s * mat[p][p];
		for (int i = 0; i < n; i++) {
			if (i != p && i != q) {
				mat[i][i] = mat[i][i];
			}
		}
	} while (max_off_diag > eps);

	// ��ȡ����ֵ
	for (int i = 0; i < n; i++) {
		eigens[i] = mat[i][i];
	}
}

void U2NITS::Space::RoeAverage(myfloat U1[4], myfloat U2[4], myfloat gamma) {
	/*
	��h1 h2
	*/
	myfloat ruvp1[4];
	U2NITS::Math::U2ruvp_host(U1, ruvp1, gamma);
	myfloat rho1 = ruvp1[0];
	myfloat u1 = ruvp1[1];
	myfloat v1 = ruvp1[2];
	myfloat p1 = ruvp1[3];
	myfloat E1 = U1[3] / U1[0];
	myfloat e1 = E1 - 0.5 * (u1 * u1 + v1 * v1);
	myfloat h1 = e1 + p1 / rho1;
	myfloat ruvp2[4];
	U2NITS::Math::U2ruvp_host(U2, ruvp2, gamma);
	myfloat rho2 = ruvp2[0];
	myfloat u2 = ruvp2[1];
	myfloat v2 = ruvp2[2];
	myfloat p2 = ruvp2[3];
	myfloat E2 = U2[3] / U2[0];
	myfloat e2 = E2 - 0.5 * (u2 * u2 + v2 * v2);
	myfloat h2 = e2 + p2 / rho2;// h = gamma*e = e+(gamma-1)*e = e+R/Cv*e = e+R*e/Cv = e+R*t = e+p/rho

	/*
	���������¼������μ�P179
	*/
	/*
	��U{j+1/2}
	Roeƽ����ʽ��u = (sr1*u1+sr2*u2)/(sr1+sr2)������sr1=sqrt(rho1), sr2=sqrt(rho2)
	�˴���	wt1 = sr1/(sr1+sr2) = rho1/(rho1+sr1*sr2) = (rho1/rho)/(rho1/rho + my_1)
	wt2 = 1-wt1 = sr2/(sr1+sr2)
	*/
	const myfloat my_half = 0.5;
	const myfloat my_1 = 1.0;
	myfloat rho = sqrt(rho1 * rho2);
	myfloat wt1 = (rho1 / rho) / (rho1 / rho + my_1);
	myfloat wt2 = my_1 - wt1;
	myfloat u = u1 * wt1 + u2 * wt2;
	myfloat v = v1 * wt1 + v2 * wt2;
	myfloat h = h1 * wt1 + h2 * wt2;
	myfloat e = h / gamma;
	myfloat p = rho * e * (gamma - 1);// ע�ⲻ��E��E=e+0.5*V2
	myfloat V2 = u * u + v * v;// �ٶ�ƽ��
	myfloat E = e + 0.5 * V2;
	myfloat U[4]{ rho,rho * u,rho * v,rho * E };// ��ά

	// ����U{j+1/2}������ֵ����������
	// // ����ֱ�Ӳ���UNITs�ķ�����
	//myfloat F[4]{ rho * u,rho * u * u + p,rho * u * v,(rho * E + p) * u };

}

void U2NITS::Space::RoeAverageFUN3D(
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
	bool if_turb1) {
	// �ο�fun3d

	myfloat wt1, wt2;// weights
	myfloat alfa_eq; // for equilibrium air
	myfloat beta_eq; // for equilibrium air
	myfloat ratio;   // weight
	myfloat dpress;  // jump in pressure
	myfloat denergy; // jump in energy
	myfloat dpbar_de;// unweighted pressure-energy derivative
	myfloat dp_denom;// denominator term Liou et. al pressure weighting
	myfloat p_res;   // residual pressure
	int ne1;      // 1st energy equation index
	int n_pjac;   // press jac elements
	const myfloat my_half = 0.5;
	const myfloat my_1 = 1.0;

	/*
	����Roeƽ����ʽ��u = (sr1*u1+sr2*u2)/(sr1+sr2)��
	����sr1=sqrt(rho1), sr2=sqrt(rho2)
	�˴���
	wt1 = sr1/(sr1+sr2) = rho1/(rho1+sr1*sr2) = (rho1/rho)/(rho1/rho + my_1)
	wt2 = 1-wt1 = sr2/(sr1+sr2)
	*/
	rho = sqrt(rho1 * rho2);
	wt1 = (rho1 / rho) / (rho1 / rho + my_1);
	wt2 = my_1 - wt1;
	u = u1 * wt1 + u2 * wt2;
	v = v1 * wt1 + v2 * wt2;
	h = h1 * wt1 + h2 * wt2;
	e = e1 * wt1 + e2 * wt2;
	frac = frac1 * wt1 + frac2 * wt2;
	if (if_turb1) {
		turb = turb1 * wt1 + turb2 * wt2;
	}

}


void U2NITS::Space::ConvectRoeCommon3d(const myfloat UL[5], const myfloat UR[5], const myfloat faceNormal[3],
	const myfloat faceArea, myfloat faceFlux[5], bool bDynamicMesh, myfloat dynamicMeshValue, myfloat gamma, myfloat rcpcv) {
	// ������������Roe����� �ο�UNITs Convect_Roe_Common
	// ��RiemannSolver����
	// �����faceFlux

	// �淨��λ����
	myfloat sav1n = faceNormal[0];
	myfloat sav2n = faceNormal[1];
	myfloat sav3n = faceNormal[2];
	myfloat sav4n = dynamicMeshValue / faceArea;//��������أ�Ŀǰ����Ҫ
	// �غ���ת������rho u v w p
	myfloat ruvwpL[5]{};
	U2NITS::Math::U2ruvwp_host_3d(UL, ruvwpL, gamma);
	myfloat ruvwpR[5]{};
	U2NITS::Math::U2ruvwp_host_3d(UR, ruvwpR, gamma);
	myfloat rL = ruvwpL[0];
	myfloat uL = ruvwpL[1];
	myfloat vL = ruvwpL[2];
	myfloat wL = ruvwpL[3];
	myfloat pL = ruvwpL[4];
	myfloat rR = ruvwpR[0];
	myfloat uR = ruvwpR[1];
	myfloat vR = ruvwpR[2];
	myfloat wR = ruvwpR[3];
	myfloat pR = ruvwpR[4];
	// �����ٶ�
	myfloat ulnormaln = (sav1n * uL + sav2n * vL + sav3n * wL + sav4n);
	myfloat urnormaln = (sav1n * uR + sav2n * vR + sav3n * wR + sav4n);
	myfloat rulnormaln = rL * ulnormaln;
	myfloat rurnormaln = rR * urnormaln;
	// �ٶ�ƽ��
	myfloat V2L = uL * uL + vL * vL + wL * wL;
	myfloat V2R = uR * uR + vR * vR + wR * wR;
	//// ��+���� H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	myfloat ga1 = gamma - 1.0;
	myfloat HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	myfloat HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// ����+���� E = e + 0.5*V2 = p/rho/(gamma-1) + 0.5*V2
	// ���������غ��� rhoE = p/(gamma-1) + 0.5*rho*V2
	myfloat rhoEL = pL / ga1 + 0.5 * rL * V2L;
	myfloat rhoER = pR / ga1 + 0.5 * rR * V2R;
	// ͨ�� Ϊʲô������H��
	auto& flux = faceFlux;
	auto& sav = faceArea;
	flux[0] = 0.5 * (rulnormaln + rurnormaln) * sav;
	flux[1] = 0.5 * (rulnormaln * uL + rurnormaln * uR + sav1n * (pL + pR)) * sav;
	flux[2] = 0.5 * (rulnormaln * vL + rurnormaln * vR + sav2n * (pL + pR)) * sav;
	flux[3] = 0.5 * (rulnormaln * wL + rurnormaln * wR + sav3n * (pL + pR)) * sav;
	flux[4] = 0.5 * (rulnormaln * HL + rurnormaln * HR - sav4n * (pL + pR)) * sav;

	// ����������
	myfloat drRoe[5]{};
	RoeDissapationTerm3d(
		gamma, UL, UR, ruvwpL, ruvwpR, faceNormal, faceArea,
		bDynamicMesh, dynamicMeshValue, 0, 0, 0, drRoe
	);

	//// ����Ӧ��ɢϵ�� funscheme���������
	//myfloat funscheme = AdaptiveFunctionCoeffient();
	//
	//myfloat droR = -(flux[0] - funscheme * drRoe[0]);
	//myfloat drxR = -(flux[1] - funscheme * drRoe[1]);
	//myfloat dryR = -(flux[2] - funscheme * drRoe[2]);
	//myfloat drzR = -(flux[3] - funscheme * drRoe[3]);
	//myfloat dreR = -(flux[4] - funscheme * drRoe[4]);

	const myfloat funscheme = 1.0;
	flux[0] -= funscheme * drRoe[0];
	flux[1] -= funscheme * drRoe[1];
	flux[2] -= funscheme * drRoe[2];
	flux[3] -= funscheme * drRoe[3];
	flux[4] -= funscheme * drRoe[4];
}

void U2NITS::Space::RoeDissapationTerm3d(
	myfloat gamma,
	const myfloat UL[5], const myfloat UR[5],
	myfloat ruvwpL[5], myfloat ruvwpR[5],
	const myfloat faceNormal[3], myfloat faceArea,
	bool bDynamicMesh, myfloat dynamicMeshValue,
	bool bEntropyFix, myfloat KEntropyFix[3], myfloat kp,
	myfloat drRoe[5]
) {
	// ����UNITs Roe_dissipation_term
	// �浥λ����
	typedef myfloat myfloat;
	myfloat ga1 = gamma - 1.0;
	myfloat sav1n = faceNormal[0];
	myfloat sav2n = faceNormal[1];
	myfloat sav3n = faceNormal[2];
	myfloat sav4n = dynamicMeshValue / faceArea;


	// ���ҵ�Ԫ������
	myfloat& rL = ruvwpL[0];
	myfloat& uL = ruvwpL[1];
	myfloat& vL = ruvwpL[2];
	myfloat& wL = ruvwpL[3];
	myfloat& pL = ruvwpL[4];
	myfloat& rR = ruvwpR[0];
	myfloat& uR = ruvwpR[1];
	myfloat& vR = ruvwpR[2];
	myfloat& wR = ruvwpR[3];
	myfloat& pR = ruvwpR[4];
	// �ٶ�ƽ��
	myfloat V2L = uL * uL + vL * vL + wL * wL;
	myfloat V2R = uR * uR + vR * vR + wR * wR;
	// ��+���� H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	myfloat HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	myfloat HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// ����
	myfloat aL = sqrt(gamma * pL / rL);
	myfloat aR = sqrt(gamma * pR / rR);
	// Roeƽ��
	myfloat rrorl = sqrt(rR / rL);
	myfloat rrorlp1 = rrorl + 1.0;
	myfloat rm = sqrt(rR * rL);
	myfloat um = (uL + uR * rrorl) / rrorlp1;
	myfloat vm = (vL + vR * rrorl) / rrorlp1;
	myfloat wm = (wL + wR * rrorl) / rrorlp1;
	myfloat Hm = (HL + HR * rrorl) / rrorlp1;
	// �����ٶ�
	myfloat vm2 = um * um + vm * vm + wm * wm;
	myfloat am2 = ga1 * (Hm - 0.5 * vm2);
	myfloat am = sqrt(am2);
	myfloat anormaln = am;
	myfloat mach2 = vm2 / am2;
	myfloat unormaln = um * sav1n + vm * sav2n + wm * sav3n; // �����ٶ�

	// ����ֵ
	myfloat eig1 = abs(unormaln);
	myfloat eig2 = abs(unormaln + anormaln);
	myfloat eig3 = abs(unormaln - anormaln);
	// ������
	const myfloat epsilon = U2NITS::Math::EPSILON;
	const myfloat& enFixK1 = KEntropyFix[0];
	const myfloat& enFixK2 = KEntropyFix[1];
	const myfloat& enFixK3 = KEntropyFix[2];
	if (bEntropyFix) {
		myfloat eig_lim1 = (abs(unormaln) + anormaln) * enFixK3 * kp;
		myfloat eig_lim2 = (abs(unormaln) + anormaln) * (enFixK1 + enFixK2 * kp);
		myfloat eig_lim3 = eig_lim2;
		// float����С����Ϊ1.4e-45
		// �μ� csappP72 �� https://zhuanlan.zhihu.com/p/656543002
		U2NITS::Space::EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
	}
	// �������̵��غ��� E=e+0.5*V2, 
	// ��+���� H = h+0.5*V2��
	// ����Ϊ e=h-p/rho��
	// ��� E=H-p/rho
	myfloat ER = HR - pR / rR;
	myfloat EL = HL - pL / rL;
	myfloat dE = ER - EL;

	// Ϊʲô������"+"����"-"��
	// ������ E=e+0.5*V2=h-p/rho+0.5*V2
	// e=h-p/rho, e=h/gamma
	myfloat Etm = Hm / gamma + ga1 / gamma * 0.5 * vm2;// ?
	// ����

	myfloat dW[5]{};
	dW[0] = rR - rL;
	dW[1] = rm * (uR - uL) + dW[0] * um;
	dW[2] = rm * (vR - vL) + dW[0] * vm;
	dW[3] = rm * (wR - wL) + dW[0] * wm;
	dW[4] = rm * (ER - EL) + dW[0] * Etm;

	myfloat astar = (eig2 + eig3) / 2.0;
	myfloat Mstar = (eig2 - eig3) / 2.0 / am;
	myfloat dunormaln = (sav1n * (uR - uL) + sav2n * (vR - vL) + sav3n * (wR - wL));
	myfloat dUroe = Mstar * dunormaln + (astar - eig1) * (pR - pL) / rm / am2;
	myfloat dProe = Mstar * (pR - pL) + (astar - eig1) * rm * dunormaln;

	// �����ɢ��
	myfloat& sav = faceArea;
	myfloat droRoe = 0.5 * (eig1 * dW[0] + dUroe * rm) * sav;
	myfloat drxRoe = 0.5 * (eig1 * dW[1] + dUroe * rm * um + dProe * sav1n) * sav;
	myfloat dryRoe = 0.5 * (eig1 * dW[2] + dUroe * rm * vm + dProe * sav2n) * sav;
	myfloat drzRoe = 0.5 * (eig1 * dW[3] + dUroe * rm * wm + dProe * sav3n) * sav;
	// unormaln - sav4n�������Ϊ�����ٶȣ�
	myfloat dreRoe = 0.5 * (eig1 * dW[4] + dUroe * rm * Hm + dProe * (unormaln - sav4n)) * sav;

	drRoe[0] = droRoe;
	drRoe[1] = drxRoe;
	drRoe[2] = dryRoe;
	drRoe[3] = drzRoe;
	drRoe[4] = dreRoe;
}

void U2NITS::Space::ConvectRoeCommon2d(const myfloat UL[4], const myfloat UR[4], const myfloat faceNormal[2], const myfloat faceArea, myfloat faceFlux[4], myfloat gamma, myfloat rcpcv) {
	// ������������Roe����� �ο�UNITs Convect_Roe_Common
	// ��RiemannSolver����
	// �����faceFlux
	myfloat nx = faceNormal[0];
	myfloat ny = faceNormal[1];
	myfloat velocity_dynaMesh = 0;//��������أ�Ŀǰ����Ҫ
	// �غ���ת������rho u v w p
	myfloat ruvpL[4]{};
	myfloat ruvpR[4]{};
	U2NITS::Math::U2ruvp_host(UL, ruvpL, gamma);
	U2NITS::Math::U2ruvp_host(UR, ruvpR, gamma);
	myfloat rL = ruvpL[0];
	myfloat uL = ruvpL[1];
	myfloat vL = ruvpL[2];
	myfloat pL = ruvpL[3];
	myfloat rR = ruvpR[0];
	myfloat uR = ruvpR[1];
	myfloat vR = ruvpR[2];
	myfloat pR = ruvpR[3];
	//// ����
	//U2NITS::Math::restrictRhoAndP(ruvpL);
	//U2NITS::Math::restrictRhoAndP(ruvpR);
	// �����ٶ�
	myfloat ulnormaln = (nx * uL + ny * vL + velocity_dynaMesh);
	myfloat urnormaln = (nx * uR + ny * vR + velocity_dynaMesh);
	myfloat rulnormaln = rL * ulnormaln;
	myfloat rurnormaln = rR * urnormaln;
	// �ٶ�ƽ��
	myfloat V2L = uL * uL + vL * vL;
	myfloat V2R = uR * uR + vR * vR;
	//// ��+���� H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	myfloat ga1 = gamma - 1.0;
	myfloat HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	myfloat HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// ����+���� E = e + 0.5*V2 = p/rho/(gamma-1) + 0.5*V2
	// ���������غ��� rhoE = p/(gamma-1) + 0.5*rho*V2
	myfloat rhoEL = pL / ga1 + 0.5 * rL * V2L;
	myfloat rhoER = pR / ga1 + 0.5 * rR * V2R;
	// ͨ�� Ϊʲô������H��
	faceFlux[0] = 0.5 * (rulnormaln + rurnormaln) * faceArea;
	faceFlux[1] = 0.5 * (rulnormaln * uL + rurnormaln * uR + nx * (pL + pR)) * faceArea;
	faceFlux[2] = 0.5 * (rulnormaln * vL + rurnormaln * vR + ny * (pL + pR)) * faceArea;
	faceFlux[3] = 0.5 * (rulnormaln * HL + rurnormaln * HR - velocity_dynaMesh * (pL + pR)) * faceArea;

	// ������ ���ڼ���̽�����ӻ�δд��Ŀǰ����������
	myfloat drRoe[4]{};
	bool entropyFix = false;
	myfloat kEntropyFix[3]{ 0.2,7.5,20.0 };// ������ϵ��(EnFix_k1 = 0.15 - 0.25, EnFix_k2 = 5.0 - 10.0, EnFix_k3 = 15.0 - 25.0)
	myfloat p_sensor = PShockWaveSensor();// ����̽�����ӣ�����ѹ���ռ���׵���������������ǿ�ȣ�ȡֵ[0,1)
	RoeDissapationTerm2d(
		gamma, ruvpL, ruvpR, faceNormal, faceArea,
		entropyFix, kEntropyFix, p_sensor, drRoe
	);


	const myfloat funscheme = 1.0;// ����Ӧ��ɢϵ��
	faceFlux[0] -= funscheme * drRoe[0];
	faceFlux[1] -= funscheme * drRoe[1];
	faceFlux[2] -= funscheme * drRoe[2];
	faceFlux[3] -= funscheme * drRoe[3];

}

void U2NITS::Space::RoeDissapationTerm2d(myfloat gamma, myfloat ruvpL[4], myfloat ruvpR[4], const myfloat faceNormal[2], myfloat faceArea, bool bEntropyFix, myfloat KEntropyFix[3], myfloat pShockWaveSensor, myfloat drRoe[4]) {
	// ����UNITs Roe_dissipation_term
	// �浥λ����

	myfloat ga1 = gamma - 1.0;
	myfloat nx = faceNormal[0];
	myfloat ny = faceNormal[1];
	myfloat meshNormalVelocity = 0.0;// ��������� �����˶��ٶ�

	// ���ҵ�Ԫ������
	myfloat rL = ruvpL[0];
	myfloat uL = ruvpL[1];
	myfloat vL = ruvpL[2];
	myfloat pL = ruvpL[3];
	myfloat rR = ruvpR[0];
	myfloat uR = ruvpR[1];
	myfloat vR = ruvpR[2];
	myfloat pR = ruvpR[3];
	// �ٶ�ƽ��
	myfloat V2L = uL * uL + vL * vL;
	myfloat V2R = uR * uR + vR * vR;
	// ��+���� H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	myfloat HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	myfloat HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// ����
	myfloat aL = sqrt(gamma * pL / rL);
	myfloat aR = sqrt(gamma * pR / rR);
	// Roeƽ�� 
	myfloat rrorl = sqrt(rR / rL);
	myfloat rrorlp1 = rrorl + 1.0;
	myfloat rm = sqrt(rR * rL);
	myfloat um = (uL + uR * rrorl) / rrorlp1;
	myfloat vm = (vL + vR * rrorl) / rrorlp1;
	myfloat Hm = (HL + HR * rrorl) / rrorlp1;
	myfloat Vm2 = um * um + vm * vm;
	myfloat am2 = ga1 * (Hm - 0.5 * Vm2);
	myfloat am = sqrt(am2);
	myfloat anormaln = am;
	myfloat mach2 = Vm2 / am2;
	// ��ת�任��ת��Ϊ����һάEuler����
	myfloat unormaln = um * nx + vm * ny; // �����ٶ�

	// ����ֵ һάEuler������3������ֵ u u+a u-a
	myfloat eig1 = abs(unormaln);
	myfloat eig2 = abs(unormaln + anormaln);
	myfloat eig3 = abs(unormaln - anormaln);
	// ������
	const myfloat epsilon = U2NITS::Math::EPSILON;
	const myfloat& enFixK1 = KEntropyFix[0];
	const myfloat& enFixK2 = KEntropyFix[1];
	const myfloat& enFixK3 = KEntropyFix[2];
	if (bEntropyFix) {
		myfloat eig_lim1 = (abs(unormaln) + anormaln) * enFixK3 * pShockWaveSensor;// kpΪ����̽�����ӣ�ȡ���൥Ԫp_sensor�ľ�ֵ
		myfloat eig_lim2 = (abs(unormaln) + anormaln) * (enFixK1 + enFixK2 * pShockWaveSensor);
		myfloat eig_lim3 = eig_lim2;

		U2NITS::Space::EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
	}
	// �������̵��غ��� E=e+0.5*V2, 
	// ��+���� H = h+0.5*V2��
	// ����Ϊ e=h-p/rho��
	// ��� E=H-p/rho
	myfloat ER = HR - pR / rR;
	myfloat EL = HL - pL / rL;
	myfloat dE = ER - EL;

	// Ϊʲô������"+"����"-"��
	// ������ E=e+0.5*V2=h-p/rho+0.5*V2
	// e=h-p/rho, e=h/gamma
	myfloat Etm = Hm / gamma + ga1 / gamma * 0.5 * Vm2;// ?
	// ����

	myfloat dW[4]{};
	dW[0] = rR - rL;
	dW[1] = rm * (uR - uL) + dW[0] * um;
	dW[2] = rm * (vR - vL) + dW[0] * vm;
	dW[3] = rm * (ER - EL) + dW[0] * Etm;

	myfloat astar = (eig2 + eig3) / 2.0;
	myfloat Mstar = (eig2 - eig3) / 2.0 / am;
	myfloat dunormaln = (nx * (uR - uL) + ny * (vR - vL));
	myfloat dUroe = Mstar * dunormaln + (astar - eig1) * (pR - pL) / rm / am2;
	myfloat dProe = Mstar * (pR - pL) + (astar - eig1) * rm * dunormaln;

	// �����ɢ�� rho rhou rhov rhoE
	myfloat& sav = faceArea;
	drRoe[0] = 0.5 * (eig1 * dW[0] + dUroe * rm) * sav;
	drRoe[1] = 0.5 * (eig1 * dW[1] + dUroe * rm * um + dProe * nx) * sav;
	drRoe[2] = 0.5 * (eig1 * dW[2] + dUroe * rm * vm + dProe * ny) * sav;
	drRoe[3] = 0.5 * (eig1 * dW[3] + dUroe * rm * Hm + dProe * (unormaln - meshNormalVelocity)) * sav;

}



void U2NITS::Space::GetRoeMatrix3d(
	myfloat gamma,
	myfloat ruvwpL[5], myfloat ruvwpR[5],
	myfloat faceVector[3], myfloat faceArea,
	bool bDynamicMesh, myfloat dynamicMeshValue,
	bool bEntropyFix, myfloat KEntropyFix[3], myfloat kp,
	myfloat roeMatrix[5][5]
) {
	// ���룺����
	// �����roeMatrix[5][5] ��������
	// ���ʣ�
	// Ϊʲô�淨������4ά����������3ά��
	// 
	typedef myfloat myfloat;// myfloat�Ǿֲ�����
	const myfloat epsilon = U2NITS::Math::EPSILON;
	myfloat& area = faceArea;
	auto& S_vector = faceVector;
	myfloat ga1 = gamma - 1;
	// �������
	area = (area < epsilon) ? epsilon : area;
	// �浥λ����������4�������Ƕ��������
	myfloat Sn[4]{};
	Sn[0] = S_vector[0] / area;
	Sn[1] = S_vector[1] / area;
	Sn[2] = S_vector[2] / area;
	Sn[3] = dynamicMeshValue / area;
	// ���ҵ�Ԫ������
	myfloat& rL = ruvwpL[0];
	myfloat& uL = ruvwpL[1];
	myfloat& vL = ruvwpL[2];
	myfloat& wL = ruvwpL[3];
	myfloat& pL = ruvwpL[4];
	myfloat& rR = ruvwpR[0];
	myfloat& uR = ruvwpR[1];
	myfloat& vR = ruvwpR[2];
	myfloat& wR = ruvwpR[3];
	myfloat& pR = ruvwpR[4];
	// �ٶ�ƽ��
	myfloat V2L = uL * uL + vL * vL + wL * wL;
	myfloat V2R = uR * uR + vR * vR + wR * wR;
	// ��+���� H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	myfloat HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	myfloat HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// ����
	myfloat aL = sqrt(gamma * pL / rL);
	myfloat aR = sqrt(gamma * pR / rR);

	// Roeƽ��
	myfloat rrorl = sqrt(rR / rL);
	myfloat rrorlp1 = rrorl + 1.0;
	myfloat rm = sqrt(rR * rL);
	myfloat um = (uL + uR * rrorl) / rrorlp1;
	myfloat vm = (vL + vR * rrorl) / rrorlp1;
	myfloat wm = (wL + wR * rrorl) / rrorlp1;
	myfloat Hm = (HL + HR * rrorl) / rrorlp1;
	// ���������
	myfloat vm2 = um * um + vm * vm + wm * wm;
	myfloat am2 = ga1 * (Hm - 0.5 * vm2);
	myfloat am = sqrt(am2);
	myfloat mach2 = vm2 / am2;
	myfloat VN = um * Sn[0] + vm * Sn[1] + wm * Sn[2]; // �����ٶ�
	myfloat machN = VN / am;

	myfloat Vec[3]{ um,vm,wm };
	myfloat VxV[3][3]{};
	// �����ԣ�U={1,1,1}'��V={1,2,3}'ʱ���˷���UxV'��Ч����ΪV'�Զ�ʶ���3x1����
	U2NITS::Matrix::mul_ixj_jxk(3, 1, 3, Vec, Vec, (myfloat*)VxV);
	myfloat SnxSn[3][3]{};
	myfloat Sn_3[3]{ Sn[0],Sn[1],Sn[2] };// ȡǰ3��
	U2NITS::Matrix::mul_ixj_jxk(3, 1, 3, Sn_3, Sn_3, (myfloat*)SnxSn);

	// ����R1L1
	myfloat R1L1[5][5]{};
	// R1L1(1,1)
	R1L1[0][0] = 1.0 - 0.5 * ga1 * mach2;// 1-0.5*(gamma-1)*Ma^2
	// R1L1(1,5)
	R1L1[0][4] = -ga1 / am2;
	// R1L1(5,1)
	R1L1[4][0] = VN * VN - 0.5 * vm2 * (1.0 + 0.5 * ga1 * mach2);
	// R1L1(5,5)
	R1L1[4][4] = -0.5 * ga1 * mach2;
	for (int i = 0; i < 3; i++) {
		// R1L1(1,2:4)
		R1L1[0][i + 1] = ga1 / am2 * Vec[i];
		// R1L1(2:4,1)
		R1L1[i + 1][0] = -(1.0 + 0.5 * ga1 * mach2) * Vec[i] + VN * Sn[i];
		// R1L1(2:4,5)
		R1L1[i + 1][4] = -ga1 / am2 * Vec[i];
		// R1L1(5,2:4)
		R1L1[4][i + 1] = (1.0 + 0.5 * ga1 * mach2) * Vec[i] - VN * Sn[i];
		// R1L1(2:4,2:4)
		for (int j = 0; j < 3; j++) {
			myfloat UNIT_i = (i == j);// ��λ����
			R1L1[i + 1][j + 1] = ga1 / am2 * VxV[i][j] + UNIT_i - SnxSn[i][j];
		}
	}

	myfloat V_Sn[3]{};
	myfloat ga1_n[3]{};
	for (int i = 0; i < 3; i++) {
		V_Sn[i] = Vec[i] - am * Sn_3[i];
		ga1_n[i] = -0.5 * ga1 / am2 * Vec[i] - 0.5 * Sn_3[i] / am;
	}
	myfloat MIDDLE1[3][3];// ��ʱ����
	U2NITS::Matrix::mul_ixj_jxk(3, 1, 3, V_Sn, ga1_n, (myfloat*)MIDDLE1);
	// ����R2L2
	myfloat R2L2[5][5]{};
	// R2L2(1,1)
	R2L2[0][0] = 0.25 * ga1 * mach2 + 0.5 * machN;
	// R2L2(1,5)
	R2L2[0][4] = 0.5 * ga1 / am2;
	// R2L2(5,1)
	R2L2[4][0] = (Hm - VN * am) * R2L2[0][0];
	// R2L2(5,5)
	R2L2[4][4] = (Hm - VN * am) * 0.5 * ga1 / am2;
	for (int i = 0; i < 3; i++) {
		// R2L2(1,2:4)
		R2L2[0][i + 1] = ga1_n[i];
		// R2L2(2:4,1)
		R2L2[i + 1][0] = V_Sn[i] * R2L2[0][0];
		// R2L2(2:4,5)
		R2L2[i + 1][4] = 0.5 * ga1 / am2 * V_Sn[i];
		// R2L2(5,2:4)
		R2L2[4][i + 1] = (Hm - VN * am) * ga1_n[i];
		// R2L2(2:4,2:4)
		for (int j = 0; j < 3; j++) {
			R2L2[i + 1][j + 1] = MIDDLE1[i][j];
		}
	}

	// ����ga1_n�������¼���
	for (int i = 0; i < 3; i++) {
		V_Sn[i] = Vec[i] + am * Sn_3[i];
		//ga1_n[i] = -0.5 * ga1 / am2 * Vec[i] - 0.5 * Sn_3[i] / am;
	}
	U2NITS::Matrix::mul_ixj_jxk(3, 1, 3, V_Sn, ga1_n, (myfloat*)MIDDLE1);
	// ����R3L3
	myfloat R3L3[5][5]{};
	// R3L3(1,1)
	R3L3[0][0] = 0.25 * ga1 * mach2 - 0.5 * machN;
	// R3L3(1,5)
	R3L3[0][4] = 0.5 * ga1 / am2;
	// R3L3(5,1)
	R3L3[4][0] = (Hm + VN * am) * R3L3[0][0];
	// R3L3(5,5)
	R3L3[4][4] = (Hm + VN * am) * 0.5 * ga1 / am2;
	for (int i = 0; i < 3; i++) {
		// R3L3(1,2:4)
		R3L3[0][i + 1] = ga1_n[i];
		// R3L3(2:4,1)
		R3L3[i + 1][0] = V_Sn[i] * R3L3[0][0];
		// R3L3(2:4,5)
		R3L3[i + 1][4] = 0.5 * ga1 / am2 * V_Sn[i];
		// R3L3(5,2:4)
		R3L3[4][i + 1] = (Hm + VN * am) * ga1_n[i];
		// R3L3(2:4,2:4)
		for (int j = 0; j < 3; j++) {
			R3L3[i + 1][j + 1] = MIDDLE1[i][j];
		}
	}

	myfloat eig1 = abs(VN + Sn[3]);
	myfloat eig2 = abs(VN + Sn[3] + am);
	myfloat eig3 = abs(VN + Sn[3] - am);
	// ������
	const myfloat& enFixK1 = KEntropyFix[0];
	const myfloat& enFixK2 = KEntropyFix[1];
	const myfloat& enFixK3 = KEntropyFix[2];
	if (bEntropyFix) {
		myfloat eig_lim1 = (abs(VN + Sn[3]) + am) * enFixK3 * kp;
		myfloat eig_lim2 = (abs(VN + Sn[3]) + am) * (enFixK1 + enFixK2 * kp);
		myfloat eig_lim3 = eig_lim2;
		// float����С����Ϊ1.4e-45
		// �μ� csappP72 �� https://zhuanlan.zhihu.com/p/656543002
		U2NITS::Space::EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
	}
	// �����������
	roeMatrix;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			roeMatrix[i][j] = 
				(eig1 * R1L1[i][j]
				+ eig3 * R2L2[i][j]
				+ eig2 * R3L3[i][j])
				* area;
		}
	}
}
