#include "RoeGPU.h"
#include "../../math/MathGPU.h"

__device__ void GPU::Space::Convection::ConvectRoeCommon2d(const myfloat UL[4], const myfloat UR[4], const myfloat faceNormal[2], const myfloat faceArea, myfloat faceFlux[4], myfloat gamma, myfloat rcpcv) {
	// 采用熵修正的Roe求解器 参考UNITs Convect_Roe_Common
	// 被RiemannSolver调用
	// 输出：faceFlux
	myfloat nx = faceNormal[0];
	myfloat ny = faceNormal[1];
	myfloat ruvpL[4]{};
	myfloat ruvpR[4]{};

	// 面法向单位向量
	myfloat velocity_dynaMesh = 0;//动网格相关，目前不需要
	// 守恒量转场变量rho u v w p
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
	// 法向速度
	myfloat ulnormaln = (nx * uL + ny * vL + velocity_dynaMesh);
	myfloat urnormaln = (nx * uR + ny * vR + velocity_dynaMesh);
	myfloat rulnormaln = rL * ulnormaln;
	myfloat rurnormaln = rR * urnormaln;
	// 速度平方
	myfloat V2L = uL * uL + vL * vL;
	myfloat V2R = uR * uR + vR * vR;
	//// 焓+动能 H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	myfloat ga1 = gamma - 1.0;
	myfloat HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	myfloat HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// 内能+动能 E = e + 0.5*V2 = p/rho/(gamma-1) + 0.5*V2
	// 能量方程守恒量 rhoE = p/(gamma-1) + 0.5*rho*V2
	myfloat rhoEL = pL / ga1 + 0.5 * rL * V2L;
	myfloat rhoER = pR / ga1 + 0.5 * rR * V2R;
	// 通量 为什么这里是H？
	faceFlux[0] = 0.5 * (rulnormaln + rurnormaln) * faceArea;
	faceFlux[1] = 0.5 * (rulnormaln * uL + rurnormaln * uR + nx * (pL + pR)) * faceArea;
	faceFlux[2] = 0.5 * (rulnormaln * vL + rurnormaln * vR + ny * (pL + pR)) * faceArea;
	faceFlux[3] = 0.5 * (rulnormaln * HL + rurnormaln * HR - velocity_dynaMesh * (pL + pR)) * faceArea;

	// 熵修正 由于激波探测因子还未写，目前禁用熵修正
	myfloat drRoe[4]{};
	bool entropyFix = false;
	myfloat kEntropyFix[3]{ 0.2,7.5,20.0 };// 熵修正系数(EnFix_k1 = 0.15 - 0.25, EnFix_k2 = 5.0 - 10.0, EnFix_k3 = 15.0 - 25.0)
	myfloat p_sensor = PShockWaveSensor();// 激波探测因子，等于压力空间二阶导。用来衡量激波强度，取值[0,1)
	RoeDissapationTerm2d(
		gamma, ruvpL, ruvpR, faceNormal, faceArea,
		entropyFix, kEntropyFix, p_sensor, drRoe
	);

	const myfloat funscheme = 1.0;// 自适应耗散系数
	faceFlux[0] -= funscheme * drRoe[0];
	faceFlux[1] -= funscheme * drRoe[1];
	faceFlux[2] -= funscheme * drRoe[2];
	faceFlux[3] -= funscheme * drRoe[3];

}

__device__ void GPU::Space::Convection::RoeDissapationTerm2d(myfloat gamma, myfloat ruvpL[4], myfloat ruvpR[4], const myfloat faceNormal[2], myfloat faceArea, bool bEntropyFix, myfloat KEntropyFix[3], myfloat pShockWaveSensor, myfloat drRoe[4]) {


	myfloat ga1 = gamma - 1.0;
	myfloat nx = faceNormal[0];
	myfloat ny = faceNormal[1];
	myfloat meshNormalVelocity = 0.0;// 动网格相关 法向运动速度

	// 左右单元场变量
	myfloat rL = ruvpL[0];
	myfloat uL = ruvpL[1];
	myfloat vL = ruvpL[2];
	myfloat pL = ruvpL[3];
	myfloat rR = ruvpR[0];
	myfloat uR = ruvpR[1];
	myfloat vR = ruvpR[2];
	myfloat pR = ruvpR[3];
	// 速度平方
	myfloat V2L = uL * uL + vL * vL;
	myfloat V2R = uR * uR + vR * vR;
	// 焓+动能 H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	myfloat HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	myfloat HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// 声速
	myfloat aL = sqrt(gamma * pL / rL);
	myfloat aR = sqrt(gamma * pR / rR);
	// Roe平均 
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
	// 旋转变换，转化为扩张一维Euler方程
	myfloat unormaln = um * nx + vm * ny; // 法向速度

	// 特征值 一维Euler方程有3个特征值 u u+a u-a
	myfloat eig1 = abs(unormaln);
	myfloat eig2 = abs(unormaln + anormaln);
	myfloat eig3 = abs(unormaln - anormaln);
	// 熵修正
	const myfloat epsilon = U2NITS::Math::EPSILON;
	const myfloat& enFixK1 = KEntropyFix[0];
	const myfloat& enFixK2 = KEntropyFix[1];
	const myfloat& enFixK3 = KEntropyFix[2];
	if (bEntropyFix) {
		myfloat eig_lim1 = (abs(unormaln) + anormaln) * enFixK3 * pShockWaveSensor;// kp为激波探测因子，取两侧单元p_sensor的均值
		myfloat eig_lim2 = (abs(unormaln) + anormaln) * (enFixK1 + enFixK2 * pShockWaveSensor);
		myfloat eig_lim3 = eig_lim2;

		EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
	}
	// 能量方程的守恒量 E=e+0.5*V2, 
	// 焓+动能 H = h+0.5*V2，
	// 又因为 e=h-p/rho，
	// 因此 E=H-p/rho
	myfloat ER = HR - pR / rR;
	myfloat EL = HL - pL / rL;
	myfloat dE = ER - EL;

	// 为什么这里是"+"不是"-"？
	// 理论上 E=e+0.5*V2=h-p/rho+0.5*V2
	// e=h-p/rho, e=h/gamma
	myfloat Etm = Hm / gamma + ga1 / gamma * 0.5 * Vm2;// ?
	// 计算

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

	// 计算耗散项 rho rhou rhov rhoE
	myfloat& sav = faceArea;
	drRoe[0] = 0.5 * (eig1 * dW[0] + dUroe * rm) * sav;
	drRoe[1] = 0.5 * (eig1 * dW[1] + dUroe * rm * um + dProe * nx) * sav;
	drRoe[2] = 0.5 * (eig1 * dW[2] + dUroe * rm * vm + dProe * ny) * sav;
	drRoe[3] = 0.5 * (eig1 * dW[3] + dUroe * rm * Hm + dProe * (unormaln - meshNormalVelocity)) * sav;


}
