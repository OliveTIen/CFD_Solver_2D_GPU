#include "RoeGPU.h"
#include "../../math/MathGPU.h"

__device__ void GPU::Space::Convection::ConvectRoeCommon2d(const REAL UL[4], const REAL UR[4], const REAL faceNormal[2], const REAL faceArea, REAL faceFlux[4], DGlobalPara& para) {
	// ������������Roe����� �ο�UNITs Convect_Roe_Common
	// ��RiemannSolver����
	// �����faceFlux
	real nx = faceNormal[0];
	real ny = faceNormal[1];
	real ruvpL[4]{};
	real ruvpR[4]{};

	// �淨��λ����
	real velocity_dynaMesh = 0;//��������أ�Ŀǰ����Ҫ
	real rcpcv = *(para.constant_R->ptr);
	real gamma = *(para.constant_gamma->ptr);
	// �غ���ת������rho u v w p
	GPU::Math::U2ruvp(UL, ruvpL, gamma);
	GPU::Math::U2ruvp(UR, ruvpR, gamma);
	real rL = ruvpL[0];
	real uL = ruvpL[1];
	real vL = ruvpL[2];
	real pL = ruvpL[3];
	real rR = ruvpR[0];
	real uR = ruvpR[1];
	real vR = ruvpR[2];
	real pR = ruvpR[3];
	// �����ٶ�
	real ulnormaln = (nx * uL + ny * vL + velocity_dynaMesh);
	real urnormaln = (nx * uR + ny * vR + velocity_dynaMesh);
	real rulnormaln = rL * ulnormaln;
	real rurnormaln = rR * urnormaln;
	// �ٶ�ƽ��
	real V2L = uL * uL + vL * vL;
	real V2R = uR * uR + vR * vR;
	//// ��+���� H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	real ga1 = gamma - 1.0;
	real HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	real HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// ����+���� E = e + 0.5*V2 = p/rho/(gamma-1) + 0.5*V2
	// ���������غ��� rhoE = p/(gamma-1) + 0.5*rho*V2
	real rhoEL = pL / ga1 + 0.5 * rL * V2L;
	real rhoER = pR / ga1 + 0.5 * rR * V2R;
	// ͨ�� Ϊʲô������H��
	faceFlux[0] = 0.5 * (rulnormaln + rurnormaln) * faceArea;
	faceFlux[1] = 0.5 * (rulnormaln * uL + rurnormaln * uR + nx * (pL + pR)) * faceArea;
	faceFlux[2] = 0.5 * (rulnormaln * vL + rurnormaln * vR + ny * (pL + pR)) * faceArea;
	faceFlux[3] = 0.5 * (rulnormaln * HL + rurnormaln * HR - velocity_dynaMesh * (pL + pR)) * faceArea;

	// ������ ���ڼ���̽�����ӻ�δд��Ŀǰ����������
	real drRoe[4]{};
	bool entropyFix = false;
	real kEntropyFix[3]{ 0.2,7.5,20.0 };// ������ϵ��(EnFix_k1 = 0.15 - 0.25, EnFix_k2 = 5.0 - 10.0, EnFix_k3 = 15.0 - 25.0)
	real p_sensor = PShockWaveSensor();// ����̽�����ӣ�����ѹ���ռ���׵���������������ǿ�ȣ�ȡֵ[0,1)
	RoeDissapationTerm2d(
		gamma, ruvpL, ruvpR, faceNormal, faceArea,
		entropyFix, kEntropyFix, p_sensor, drRoe
	);

	const real funscheme = 1.0;// ����Ӧ��ɢϵ��
	faceFlux[0] -= funscheme * drRoe[0];
	faceFlux[1] -= funscheme * drRoe[1];
	faceFlux[2] -= funscheme * drRoe[2];
	faceFlux[3] -= funscheme * drRoe[3];

}

__device__ void GPU::Space::Convection::RoeDissapationTerm2d(REAL gamma, REAL ruvpL[4], REAL ruvpR[4], const REAL faceNormal[2], REAL faceArea, bool bEntropyFix, REAL KEntropyFix[3], REAL pShockWaveSensor, REAL drRoe[4]) {


	real ga1 = gamma - 1.0;
	real nx = faceNormal[0];
	real ny = faceNormal[1];
	real meshNormalVelocity = 0.0;// ��������� �����˶��ٶ�

	// ���ҵ�Ԫ������
	real rL = ruvpL[0];
	real uL = ruvpL[1];
	real vL = ruvpL[2];
	real pL = ruvpL[3];
	real rR = ruvpR[0];
	real uR = ruvpR[1];
	real vR = ruvpR[2];
	real pR = ruvpR[3];
	// �ٶ�ƽ��
	real V2L = uL * uL + vL * vL;
	real V2R = uR * uR + vR * vR;
	// ��+���� H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	real HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	real HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// ����
	real aL = sqrt(gamma * pL / rL);
	real aR = sqrt(gamma * pR / rR);
	// Roeƽ�� 
	real rrorl = sqrt(rR / rL);
	real rrorlp1 = rrorl + 1.0;
	real rm = sqrt(rR * rL);
	real um = (uL + uR * rrorl) / rrorlp1;
	real vm = (vL + vR * rrorl) / rrorlp1;
	real Hm = (HL + HR * rrorl) / rrorlp1;
	real Vm2 = um * um + vm * vm;
	real am2 = ga1 * (Hm - 0.5 * Vm2);
	real am = sqrt(am2);
	real anormaln = am;
	real mach2 = Vm2 / am2;
	// ��ת�任��ת��Ϊ����һάEuler����
	real unormaln = um * nx + vm * ny; // �����ٶ�

	// ����ֵ һάEuler������3������ֵ u u+a u-a
	real eig1 = abs(unormaln);
	real eig2 = abs(unormaln + anormaln);
	real eig3 = abs(unormaln - anormaln);
	// ������
	const real epsilon = U2NITS::Math::EPSILON;
	const real& enFixK1 = KEntropyFix[0];
	const real& enFixK2 = KEntropyFix[1];
	const real& enFixK3 = KEntropyFix[2];
	if (bEntropyFix) {
		real eig_lim1 = (abs(unormaln) + anormaln) * enFixK3 * pShockWaveSensor;// kpΪ����̽�����ӣ�ȡ���൥Ԫp_sensor�ľ�ֵ
		real eig_lim2 = (abs(unormaln) + anormaln) * (enFixK1 + enFixK2 * pShockWaveSensor);
		real eig_lim3 = eig_lim2;

		EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
	}
	// �������̵��غ��� E=e+0.5*V2, 
	// ��+���� H = h+0.5*V2��
	// ����Ϊ e=h-p/rho��
	// ��� E=H-p/rho
	real ER = HR - pR / rR;
	real EL = HL - pL / rL;
	real dE = ER - EL;

	// Ϊʲô������"+"����"-"��
	// ������ E=e+0.5*V2=h-p/rho+0.5*V2
	// e=h-p/rho, e=h/gamma
	real Etm = Hm / gamma + ga1 / gamma * 0.5 * Vm2;// ?
	// ����

	real dW[4]{};
	dW[0] = rR - rL;
	dW[1] = rm * (uR - uL) + dW[0] * um;
	dW[2] = rm * (vR - vL) + dW[0] * vm;
	dW[3] = rm * (ER - EL) + dW[0] * Etm;

	real astar = (eig2 + eig3) / 2.0;
	real Mstar = (eig2 - eig3) / 2.0 / am;
	real dunormaln = (nx * (uR - uL) + ny * (vR - vL));
	real dUroe = Mstar * dunormaln + (astar - eig1) * (pR - pL) / rm / am2;
	real dProe = Mstar * (pR - pL) + (astar - eig1) * rm * dunormaln;

	// �����ɢ�� rho rhou rhov rhoE
	real& sav = faceArea;
	drRoe[0] = 0.5 * (eig1 * dW[0] + dUroe * rm) * sav;
	drRoe[1] = 0.5 * (eig1 * dW[1] + dUroe * rm * um + dProe * nx) * sav;
	drRoe[2] = 0.5 * (eig1 * dW[2] + dUroe * rm * vm + dProe * ny) * sav;
	drRoe[3] = 0.5 * (eig1 * dW[3] + dUroe * rm * Hm + dProe * (unormaln - meshNormalVelocity)) * sav;


}
