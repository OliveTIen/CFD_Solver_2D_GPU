#include "RoeGPU.h"
#include "../../math/MathGPU.h"

__device__ void GPU::Space::Convection::ConvectRoeCommon2d(const myfloat UL[4], const myfloat UR[4], const myfloat faceNormal[2], const myfloat faceArea, myfloat faceFlux[4], myfloat gamma, myfloat rcpcv) {
	// ������������Roe����� �ο�UNITs Convect_Roe_Common
	// ��RiemannSolver����
	// �����faceFlux
	myfloat nx = faceNormal[0];
	myfloat ny = faceNormal[1];
	myfloat ruvpL[4]{};
	myfloat ruvpR[4]{};

	// �淨��λ����
	myfloat velocity_dynaMesh = 0;//��������أ�Ŀǰ����Ҫ
	// �غ���ת������rho u v w p
	GPU::Math::U2ruvp_device(UL, ruvpL, gamma);
	GPU::Math::U2ruvp_device(UR, ruvpR, gamma);
	myfloat rL = ruvpL[0];
	myfloat uL = ruvpL[1];
	myfloat vL = ruvpL[2];
	myfloat pL = ruvpL[3];
	myfloat rR = ruvpR[0];
	myfloat uR = ruvpR[1];
	myfloat vR = ruvpR[2];
	myfloat pR = ruvpR[3];
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

__device__ void GPU::Space::Convection::RoeDissapationTerm2d(myfloat gamma, myfloat ruvpL[4], myfloat ruvpR[4], const myfloat faceNormal[2], myfloat faceArea, bool bEntropyFix, myfloat KEntropyFix[3], myfloat pShockWaveSensor, myfloat drRoe[4]) {


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

		EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
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
