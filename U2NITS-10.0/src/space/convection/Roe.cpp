#include "Roe.h"
#include "../../math/PhysicalConvertKernel.h"
#include "../../math/GPUMatrixKernel.h"
#include "../../global/GlobalPara.h"
#include <cmath>

void U2NITS::Space::EigenValueAndVector4x4(REAL mat[4][4], REAL eigens[4], REAL R[4][4]) {
	// 计算4x4矩阵的特征值和特征向量
	// 建议直接参照UNITs的方法求
}

void U2NITS::Space::JacobiMethod(REAL mat[4][4], REAL eigens[4], REAL R[4][4]) {
	// 将R初始化为单位矩阵
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
	const REAL eps = 1.0e-10; // 收敛条件
	REAL max_off_diag; // 最大非对角元素
	REAL theta, t, c, s; // 旋转角度和旋转矩阵
	int p, q;
	do {
		max_off_diag = 0.0;
		// 寻找最大非对角元素
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
			break; // 已经足够接近对角化
		}

		// 计算旋转角度
		theta = 0.5 * std::atan2(2 * mat[p][q], mat[p][p] - mat[q][q]);// atan2返回方位角
		t = mat[p][q];
		c = std::cos(theta);
		s = std::sin(theta);

		// 修改矩阵和右特征向量
		mat[p][q] = 0.0;
		for (int i = 0; i < n; i++) {
			if (i != p && i != q) {
				REAL temp = mat[i][p];
				mat[i][p] = temp * c - mat[i][q] * s;
				mat[p][i] = mat[i][p];
				mat[i][q] = mat[i][q] * c + temp * s;
				mat[q][i] = mat[i][q];
			}
			// 更新右特征向量
			REAL new_R_ip = R[i][p] * c - R[i][q] * s;
			REAL new_R_iq = R[i][q] * c + R[i][p] * s;
			R[i][p] = new_R_ip;
			R[i][q] = new_R_iq;
		}
		mat[p][p] = mat[p][p] * c * c - 2 * t * c * s + mat[q][q] * s * s;
		mat[q][q] = mat[q][q] * c * c + 2 * t * c * s + mat[p][p] * s * s;
		mat[p][p] = mat[q][q]; // 调整无限小元素
		mat[q][q] = c * c * mat[q][q] + 2 * t * c * s + s * s * mat[p][p];
		for (int i = 0; i < n; i++) {
			if (i != p && i != q) {
				mat[i][i] = mat[i][i];
			}
		}
	} while (max_off_diag > eps);

	// 提取特征值
	for (int i = 0; i < n; i++) {
		eigens[i] = mat[i][i];
	}
}

void U2NITS::Space::RoeAverage(REAL U1[4], REAL U2[4], REAL gamma) {
	/*
	求h1 h2 
	*/
	REAL ruvp1[4];
	U2NITS::Math::U2ruvp_host(U1, ruvp1, gamma);
	REAL rho1 = ruvp1[0];
	REAL u1 = ruvp1[1];
	REAL v1 = ruvp1[2];
	REAL p1 = ruvp1[3];
	REAL E1 = U1[3] / U1[0];
	REAL e1 = E1 - 0.5 * (u1 * u1 + v1 * v1);
	REAL h1 = e1 + p1 / rho1;
	REAL ruvp2[4];
	U2NITS::Math::U2ruvp_host(U2, ruvp2, gamma);
	REAL rho2 = ruvp2[0];
	REAL u2 = ruvp2[1];
	REAL v2 = ruvp2[2];
	REAL p2 = ruvp2[3];
	REAL E2 = U2[3] / U2[0];
	REAL e2 = E2 - 0.5 * (u2 * u2 + v2 * v2);
	REAL h2 = e2 + p2 / rho2;// h = gamma*e = e+(gamma-1)*e = e+R/Cv*e = e+R*e/Cv = e+R*t = e+p/rho

	/*
	根据任玉新计算流课件P179
	*/
	/*
	求U{j+1/2}
	Roe平均公式是u = (sr1*u1+sr2*u2)/(sr1+sr2)，其中sr1=sqrt(rho1), sr2=sqrt(rho2)
	此处令	wt1 = sr1/(sr1+sr2) = rho1/(rho1+sr1*sr2) = (rho1/rho)/(rho1/rho + my_1)
	wt2 = 1-wt1 = sr2/(sr1+sr2)
	*/
	const REAL my_half = 0.5;
	const REAL my_1 = 1.0;
	REAL rho = sqrt(rho1 * rho2);
	REAL wt1 = (rho1 / rho) / (rho1 / rho + my_1);
	REAL wt2 = my_1 - wt1;
	REAL u = u1 * wt1 + u2 * wt2;
	REAL v = v1 * wt1 + v2 * wt2;
	REAL h = h1 * wt1 + h2 * wt2;
	REAL e = h / gamma;
	REAL p = rho * e * (gamma - 1);// 注意不是E。E=e+0.5*V2
	REAL V2 = u * u + v * v;// 速度平方
	REAL E = e + 0.5 * V2;
	REAL U[4]{ rho,rho * u,rho * v,rho * E };// 二维

	// 计算U{j+1/2}的特征值和特征向量
	// // 建议直接参照UNITs的方法求
	//REAL F[4]{ rho * u,rho * u * u + p,rho * u * v,(rho * E + p) * u };

}

void U2NITS::Space::RoeAverageFUN3D(
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
	bool if_turb1){
// 参考fun3d

	REAL wt1, wt2;// weights
	REAL alfa_eq; // for equilibrium air
	REAL beta_eq; // for equilibrium air
	REAL ratio;   // weight
	REAL dpress;  // jump in pressure
	REAL denergy; // jump in energy
	REAL dpbar_de;// unweighted pressure-energy derivative
	REAL dp_denom;// denominator term Liou et. al pressure weighting
	REAL p_res;   // residual pressure
	int ne1;      // 1st energy equation index
	int n_pjac;   // press jac elements
	const REAL my_half = 0.5;
	const REAL my_1 = 1.0;

	/*
	本来Roe平均公式是u = (sr1*u1+sr2*u2)/(sr1+sr2)，
	其中sr1=sqrt(rho1), sr2=sqrt(rho2)
	此处令
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


void U2NITS::Space::ConvectRoeCommon3d(const REAL UL[5], const REAL UR[5], const REAL faceNormal[3],
	const REAL faceArea, REAL faceFlux[5], bool bDynamicMesh, REAL dynamicMeshValue, REAL gamma) {
	// 采用熵修正的Roe求解器 参考UNITs Convect_Roe_Common
	// 被RiemannSolver调用
	// 输出：faceFlux

	// 面法向单位向量
	real sav1n = faceNormal[0];
	real sav2n = faceNormal[1];
	real sav3n = faceNormal[2];
	real sav4n = dynamicMeshValue / faceArea;//动网格相关，目前不需要
	real rcpcv = GlobalPara::constant::R;
	// 守恒量转场变量rho u v w p
	real ruvwpL[5]{};
	U2NITS::Math::U2ruvwp_host(UL, ruvwpL, gamma);
	real ruvwpR[5]{};
	U2NITS::Math::U2ruvwp_host(UR, ruvwpR, gamma);
	real rL = ruvwpL[0];
	real uL = ruvwpL[1];
	real vL = ruvwpL[2];
	real wL = ruvwpL[3];
	real pL = ruvwpL[4];
	real rR = ruvwpR[0];
	real uR = ruvwpR[1];
	real vR = ruvwpR[2];
	real wR = ruvwpR[3];
	real pR = ruvwpR[4];
	// 法向速度
	real ulnormaln = (sav1n * uL + sav2n * vL + sav3n * wL + sav4n);
	real urnormaln = (sav1n * uR + sav2n * vR + sav3n * wR + sav4n);
	real rulnormaln = rL * ulnormaln;
	real rurnormaln = rR * urnormaln;
	// 速度平方
	real V2L = uL * uL + vL * vL + wL * wL;
	real V2R = uR * uR + vR * vR + wR * wR;
	//// 焓+动能 H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	real ga1 = gamma - 1.0;
	real HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	real HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// 内能+动能 E = e + 0.5*V2 = p/rho/(gamma-1) + 0.5*V2
	// 能量方程守恒量 rhoE = p/(gamma-1) + 0.5*rho*V2
	real rhoEL = pL / ga1 + 0.5 * rL * V2L;
	real rhoER = pR / ga1 + 0.5 * rR * V2R;
	// 通量 为什么这里是H？
	auto& flux = faceFlux;
	auto& sav = faceArea;
	flux[0] = 0.5 * (rulnormaln + rurnormaln) * sav;
	flux[1] = 0.5 * (rulnormaln * uL + rurnormaln * uR + sav1n * (pL + pR)) * sav;
	flux[2] = 0.5 * (rulnormaln * vL + rurnormaln * vR + sav2n * (pL + pR)) * sav;
	flux[3] = 0.5 * (rulnormaln * wL + rurnormaln * wR + sav3n * (pL + pR)) * sav;
	flux[4] = 0.5 * (rulnormaln * HL + rurnormaln * HR - sav4n * (pL + pR)) * sav;

	// 不含熵修正
	real drRoe[5]{};
	RoeDissapationTerm3d(
		gamma, UL, UR, ruvwpL, ruvwpR, faceNormal, faceArea,
		bDynamicMesh, dynamicMeshValue, 0, 0, 0, drRoe
	);

	//// 自适应耗散系数 funscheme在这里添加
	//real funscheme = AdaptiveFunctionCoeffient();
	//
	//real droR = -(flux[0] - funscheme * drRoe[0]);
	//real drxR = -(flux[1] - funscheme * drRoe[1]);
	//real dryR = -(flux[2] - funscheme * drRoe[2]);
	//real drzR = -(flux[3] - funscheme * drRoe[3]);
	//real dreR = -(flux[4] - funscheme * drRoe[4]);

	const real funscheme = 1.0;
	flux[0] -= funscheme * drRoe[0];
	flux[1] -= funscheme * drRoe[1];
	flux[2] -= funscheme * drRoe[2];
	flux[3] -= funscheme * drRoe[3];
	flux[4] -= funscheme * drRoe[4];
}

void U2NITS::Space::RoeDissapationTerm3d(
	REAL gamma,
	const REAL UL[5], const REAL UR[5],
	REAL ruvwpL[5], REAL ruvwpR[5],
	const REAL faceNormal[3], REAL faceArea,
	bool bDynamicMesh, REAL dynamicMeshValue,
	bool bEntropyFix, REAL KEntropyFix[3], REAL kp,

	REAL drRoe[5]
) {
	// 参照UNITs Roe_dissipation_term
	// 面单位法向
	typedef REAL real;
	real ga1 = gamma - 1.0;
	real sav1n = faceNormal[0];
	real sav2n = faceNormal[1];
	real sav3n = faceNormal[2];
	real sav4n = dynamicMeshValue / faceArea;
	

	// 左右单元场变量
	real& rL = ruvwpL[0];
	real& uL = ruvwpL[1];
	real& vL = ruvwpL[2];
	real& wL = ruvwpL[3];
	real& pL = ruvwpL[4];
	real& rR = ruvwpR[0];
	real& uR = ruvwpR[1];
	real& vR = ruvwpR[2];
	real& wR = ruvwpR[3];
	real& pR = ruvwpR[4];
	// 速度平方
	real V2L = uL * uL + vL * vL + wL * wL;
	real V2R = uR * uR + vR * vR + wR * wR;
	// 焓+动能 H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	real HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	real HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// 声速
	real aL = sqrt(gamma * pL / rL);
	real aR = sqrt(gamma * pR / rR);
	// Roe平均
	real rrorl = sqrt(rR / rL);
	real rrorlp1 = rrorl + 1.0;
	real rm = sqrt(rR * rL);
	real um = (uL + uR * rrorl) / rrorlp1;
	real vm = (vL + vR * rrorl) / rrorlp1;
	real wm = (wL + wR * rrorl) / rrorlp1;
	real Hm = (HL + HR * rrorl) / rrorlp1;
	// 求法向速度
	real vm2 = um * um + vm * vm + wm * wm;
	real am2 = ga1 * (Hm - 0.5 * vm2);
	real am = sqrt(am2);
	real anormaln = am;
	real mach2 = vm2 / am2;
	real unormaln = um * sav1n + vm * sav2n + wm * sav3n; // 法向速度

	// 特征值
	real eig1 = abs(unormaln);
	real eig2 = abs(unormaln + anormaln);
	real eig3 = abs(unormaln - anormaln);
	// 熵修正
	const real epsilon = GPU::Matrix::EPSILON;
	const real& enFixK1 = KEntropyFix[0];
	const real& enFixK2 = KEntropyFix[1];
	const real& enFixK3 = KEntropyFix[2];
	if (bEntropyFix) {
		real eig_lim1 = (abs(unormaln) + anormaln) * enFixK3 * kp;
		real eig_lim2 = (abs(unormaln) + anormaln) * (enFixK1 + enFixK2 * kp);
		real eig_lim3 = eig_lim2;
		// float的最小正数为1.4e-45
		// 参见 csappP72 或 https://zhuanlan.zhihu.com/p/656543002
		U2NITS::Space::EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
	}
	// 能量方程的守恒量 E=e+0.5*V2, 
	// 焓+动能 H = h+0.5*V2，
	// 又因为 e=h-p/rho，
	// 因此 E=H-p/rho
	real ER = HR - pR / rR;
	real EL = HL - pL / rL;
	real dE = ER - EL;

	// 为什么这里是"+"不是"-"？
	// 理论上 E=e+0.5*V2=h-p/rho+0.5*V2
	// e=h-p/rho, e=h/gamma
	real Etm = Hm / gamma + ga1 / gamma * 0.5 * vm2;// ?
	// 计算

	real dW[5]{};
	dW[0] = rR - rL;
	dW[1] = rm * (uR - uL) + dW[0] * um;
	dW[2] = rm * (vR - vL) + dW[0] * vm;
	dW[3] = rm * (wR - wL) + dW[0] * wm;
	dW[4] = rm * (ER - EL) + dW[0] * Etm;

	real astar = (eig2 + eig3) / 2.0;
	real Mstar = (eig2 - eig3) / 2.0 / am;
	real dunormaln = (sav1n * (uR - uL) + sav2n * (vR - vL) + sav3n * (wR - wL));
	real dUroe = Mstar * dunormaln + (astar - eig1) * (pR - pL) / rm / am2;
	real dProe = Mstar * (pR - pL) + (astar - eig1) * rm * dunormaln;

	// 计算耗散项
	real& sav = faceArea;
	real droRoe = 0.5 * (eig1 * dW[0] + dUroe * rm) * sav;
	real drxRoe = 0.5 * (eig1 * dW[1] + dUroe * rm * um + dProe * sav1n) * sav;
	real dryRoe = 0.5 * (eig1 * dW[2] + dUroe * rm * vm + dProe * sav2n) * sav;
	real drzRoe = 0.5 * (eig1 * dW[3] + dUroe * rm * wm + dProe * sav3n) * sav;
	// unormaln - sav4n可以理解为法向速度，
	real dreRoe = 0.5 * (eig1 * dW[4] + dUroe * rm * Hm + dProe * (unormaln - sav4n)) * sav;

	drRoe[0] = droRoe;
	drRoe[1] = drxRoe;
	drRoe[2] = dryRoe;
	drRoe[3] = drzRoe;
	drRoe[4] = dreRoe;
}

void U2NITS::Space::ConvectRoeCommon2d(const REAL UL[4], const REAL UR[4], const REAL faceNormal[2], const REAL faceArea, REAL faceFlux[4], REAL gamma) {
	// 采用熵修正的Roe求解器 参考UNITs Convect_Roe_Common
	// 被RiemannSolver调用
	// 输出：faceFlux
	real nx = faceNormal[0];
	real ny = faceNormal[1];
	//real* ruvpL;
	//real* ruvpR;
	//ruvpL = new real[4];
	//ruvpR = new real[4];
	real ruvpL[4]{};
	real ruvpR[4]{};

	// 面法向单位向量
	real velocity_dynaMesh = 0;//动网格相关，目前不需要
	real rcpcv = GlobalPara::constant::R;
	// 守恒量转场变量rho u v w p
	U2NITS::Math::U2ruvp_host(UL, ruvpL, gamma);
	U2NITS::Math::U2ruvp_host(UR, ruvpR, gamma);
	real rL = ruvpL[0];
	real uL = ruvpL[1];
	real vL = ruvpL[2];
	real pL = ruvpL[3];
	real rR = ruvpR[0];
	real uR = ruvpR[1];
	real vR = ruvpR[2];
	real pR = ruvpR[3];
	// 法向速度
	real ulnormaln = (nx * uL + ny * vL + velocity_dynaMesh);
	real urnormaln = (nx * uR + ny * vR + velocity_dynaMesh);
	real rulnormaln = rL * ulnormaln;
	real rurnormaln = rR * urnormaln;
	// 速度平方
	real V2L = uL * uL + vL * vL;
	real V2R = uR * uR + vR * vR;
	//// 焓+动能 H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	real ga1 = gamma - 1.0;
	real HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	real HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// 内能+动能 E = e + 0.5*V2 = p/rho/(gamma-1) + 0.5*V2
	// 能量方程守恒量 rhoE = p/(gamma-1) + 0.5*rho*V2
	real rhoEL = pL / ga1 + 0.5 * rL * V2L;
	real rhoER = pR / ga1 + 0.5 * rR * V2R;
	// 通量 为什么这里是H？
	faceFlux[0] = 0.5 * (rulnormaln + rurnormaln) * faceArea;
	faceFlux[1] = 0.5 * (rulnormaln * uL + rurnormaln * uR + nx * (pL + pR)) * faceArea;
	faceFlux[2] = 0.5 * (rulnormaln * vL + rurnormaln * vR + ny * (pL + pR)) * faceArea;
	faceFlux[3] = 0.5 * (rulnormaln * HL + rurnormaln * HR - velocity_dynaMesh * (pL + pR)) * faceArea;

	// 熵修正 由于激波探测因子还未写，目前禁用熵修正
	real drRoe[4]{};
	bool entropyFix = false;
	real kEntropyFix[3]{ 0.2,7.5,20.0 };// 熵修正系数(EnFix_k1 = 0.15 - 0.25, EnFix_k2 = 5.0 - 10.0, EnFix_k3 = 15.0 - 25.0)
	real p_sensor = PShockWaveSensor();// 激波探测因子，等于压力空间二阶导。用来衡量激波强度，取值[0,1)
	RoeDissapationTerm2d(
		gamma,  ruvpL, ruvpR, faceNormal, faceArea,
		entropyFix, kEntropyFix, p_sensor, drRoe
	);


	const real funscheme = 1.0;// 自适应耗散系数
	faceFlux[0] -= funscheme * drRoe[0];
	faceFlux[1] -= funscheme * drRoe[1];
	faceFlux[2] -= funscheme * drRoe[2];
	faceFlux[3] -= funscheme * drRoe[3];
	// 在debug模式下，函数结束时会出现 Runtime check failure #2 the stack around ... was corrupted 问题
	// 但release模式下没问题
	// 将ruvpL和ruvpR设置为堆指针后，出现死循环，检查发现停留在delete[] ruvpL上
	// 后缀无论是cu还是cpp都一样
	// 而且“复制详细信息”时显示open clipboard failed
	// 已解决：原来是U2ruvwp_host函数没有修改
}

void U2NITS::Space::RoeDissapationTerm2d(REAL gamma, REAL ruvpL[4], REAL ruvpR[4], const REAL faceNormal[2], REAL faceArea, bool bEntropyFix, REAL KEntropyFix[3], REAL pShockWaveSensor, REAL drRoe[4]) {
	// 参照UNITs Roe_dissipation_term
	// 面单位法向

	real ga1 = gamma - 1.0;
	real nx = faceNormal[0];
	real ny = faceNormal[1];
	real meshNormalVelocity = 0.0;// 动网格相关 法向运动速度

	// 左右单元场变量
	real rL = ruvpL[0];
	real uL = ruvpL[1];
	real vL = ruvpL[2];
	real pL = ruvpL[3];
	real rR = ruvpR[0];
	real uR = ruvpR[1];
	real vR = ruvpR[2];
	real pR = ruvpR[3];
	// 速度平方
	real V2L = uL * uL + vL * vL;
	real V2R = uR * uR + vR * vR;
	// 焓+动能 H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	real HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	real HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// 声速
	real aL = sqrt(gamma * pL / rL);
	real aR = sqrt(gamma * pR / rR);
	// Roe平均 
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
	// 旋转变换，转化为扩张一维Euler方程
	real unormaln = um * nx + vm * ny; // 法向速度

	// 特征值 一维Euler方程有3个特征值 u u+a u-a
	real eig1 = abs(unormaln);
	real eig2 = abs(unormaln + anormaln);
	real eig3 = abs(unormaln - anormaln);
	// 熵修正
	const real epsilon = GPU::Matrix::EPSILON;
	const real& enFixK1 = KEntropyFix[0];
	const real& enFixK2 = KEntropyFix[1];
	const real& enFixK3 = KEntropyFix[2];
	if (bEntropyFix) {
		real eig_lim1 = (abs(unormaln) + anormaln) * enFixK3 * pShockWaveSensor;// kp为激波探测因子，取两侧单元p_sensor的均值
		real eig_lim2 = (abs(unormaln) + anormaln) * (enFixK1 + enFixK2 * pShockWaveSensor);
		real eig_lim3 = eig_lim2;
		
		U2NITS::Space::EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
	}
	// 能量方程的守恒量 E=e+0.5*V2, 
	// 焓+动能 H = h+0.5*V2，
	// 又因为 e=h-p/rho，
	// 因此 E=H-p/rho
	real ER = HR - pR / rR;
	real EL = HL - pL / rL;
	real dE = ER - EL;

	// 为什么这里是"+"不是"-"？
	// 理论上 E=e+0.5*V2=h-p/rho+0.5*V2
	// e=h-p/rho, e=h/gamma
	real Etm = Hm / gamma + ga1 / gamma * 0.5 * Vm2;// ?
	// 计算

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

	// 计算耗散项 rho rhou rhov rhoE
	real& sav = faceArea;
	drRoe[0] = 0.5 * (eig1 * dW[0] + dUroe * rm) * sav;
	drRoe[1] = 0.5 * (eig1 * dW[1] + dUroe * rm * um + dProe * nx) * sav;
	drRoe[2] = 0.5 * (eig1 * dW[2] + dUroe * rm * vm + dProe * ny) * sav;
	drRoe[3] = 0.5 * (eig1 * dW[3] + dUroe * rm * Hm + dProe * (unormaln - meshNormalVelocity)) * sav;

}



void U2NITS::Space::GetRoeMatrix3d(
	REAL gamma,
	REAL ruvwpL[5], REAL ruvwpR[5],
	REAL faceVector[3], REAL faceArea,
	bool bDynamicMesh, REAL dynamicMeshValue,
	bool bEntropyFix, REAL KEntropyFix[3], REAL kp,
	REAL roeMatrix[5][5]
) {
	// 输入：左右
	// 输出：roeMatrix[5][5] 特征矩阵
	// 疑问：
	// 为什么面法向量是4维向量而不是3维？
	// 
	typedef REAL real;// real是局部变量
	const real epsilon = GPU::Matrix::EPSILON;
	real& area = faceArea;
	auto& S_vector = faceVector;
	real ga1 = gamma - 1;
	// 修正面积
	area = (area < epsilon) ? epsilon : area;
	// 面单位法向量。第4个分量是动网格参数
	real Sn[4]{};
	Sn[0] = S_vector[0] / area;
	Sn[1] = S_vector[1] / area;
	Sn[2] = S_vector[2] / area;
	Sn[3] = dynamicMeshValue / area;
	// 左右单元场变量
	real& rL = ruvwpL[0];
	real& uL = ruvwpL[1];
	real& vL = ruvwpL[2];
	real& wL = ruvwpL[3];
	real& pL = ruvwpL[4];
	real& rR = ruvwpR[0];
	real& uR = ruvwpR[1];
	real& vR = ruvwpR[2];
	real& wR = ruvwpR[3];
	real& pR = ruvwpR[4];
	// 速度平方
	real V2L = uL * uL + vL * vL + wL * wL;
	real V2R = uR * uR + vR * vR + wR * wR;
	// 焓+动能 H = h + 0.5*V2 = gamma*e + 0.5*V2 =  gamma*p/rho/(gamma-1) + 0.5*V2
	real HL = pL * gamma / (rL * ga1) + 0.5 * V2L;
	real HR = pR * gamma / (rR * ga1) + 0.5 * V2R;
	// 声速
	real aL = sqrt(gamma * pL / rL);
	real aR = sqrt(gamma * pR / rR);
	
	// Roe平均
	real rrorl = sqrt(rR / rL);
	real rrorlp1 = rrorl + 1.0;
	real rm = sqrt(rR * rL);
	real um = (uL + uR * rrorl) / rrorlp1;
	real vm = (vL + vR * rrorl) / rrorlp1;
	real wm = (wL + wR * rrorl) / rrorlp1;
	real Hm = (HL + HR * rrorl) / rrorlp1;
	// 求法向马赫数
	real vm2 = um * um + vm * vm + wm * wm;
	real am2 = ga1 * (Hm - 0.5 * vm2);
	real am = sqrt(am2);
	real mach2 = vm2 / am2;
	real VN = um * Sn[0] + vm * Sn[1] + wm * Sn[2]; // 法向速度
	real machN = VN / am;

	real Vec[3]{ um,vm,wm };
	real VxV[3][3]{};
	// 经测试，U={1,1,1}'和V={1,2,3}'时，乘法对UxV'有效，因为V'自动识别成3x1矩阵
	U2NITS::Matrix::mul_ixj_jxk(3, 1, 3, Vec, Vec, (real*)VxV);
	real SnxSn[3][3]{};
	real Sn_3[3]{ Sn[0],Sn[1],Sn[2] };// 取前3项
	U2NITS::Matrix::mul_ixj_jxk(3, 1, 3, Sn_3, Sn_3, (real*)SnxSn);

	// 计算R1L1
	real R1L1[5][5]{};
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
			real UNIT_i = (i == j);// 单位矩阵
			R1L1[i + 1][j + 1] = ga1 / am2 * VxV[i][j] + UNIT_i - SnxSn[i][j];
		}
	}

	real V_Sn[3]{};
	real ga1_n[3]{};
	for (int i = 0; i < 3; i++) {
		V_Sn[i] = Vec[i] - am * Sn_3[i];
		ga1_n[i] = -0.5 * ga1 / am2 * Vec[i] - 0.5 * Sn_3[i] / am;
	}
	real MIDDLE1[3][3];// 临时变量
	U2NITS::Matrix::mul_ixj_jxk(3, 1, 3, V_Sn, ga1_n, (real*)MIDDLE1);
	// 计算R2L2
	real R2L2[5][5]{};
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

	// 这里ga1_n无需重新计算
	for (int i = 0; i < 3; i++) {
		V_Sn[i] = Vec[i] + am * Sn_3[i];
		//ga1_n[i] = -0.5 * ga1 / am2 * Vec[i] - 0.5 * Sn_3[i] / am;
	}
	U2NITS::Matrix::mul_ixj_jxk(3, 1, 3, V_Sn, ga1_n, (real*)MIDDLE1);
	// 计算R3L3
	real R3L3[5][5]{};
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

	real eig1 = abs(VN + Sn[3]);
	real eig2 = abs(VN + Sn[3] + am);
	real eig3 = abs(VN + Sn[3] - am);
	// 熵修正
	const real& enFixK1 = KEntropyFix[0];
	const real& enFixK2 = KEntropyFix[1];
	const real& enFixK3 = KEntropyFix[2];
	if (bEntropyFix) {
		real eig_lim1 = (abs(VN + Sn[3]) + am) * enFixK3 * kp;
		real eig_lim2 = (abs(VN + Sn[3]) + am) * (enFixK1 + enFixK2 * kp);
		real eig_lim3 = eig_lim2;
		// float的最小正数为1.4e-45
		// 参见 csappP72 或 https://zhuanlan.zhihu.com/p/656543002
		U2NITS::Space::EigenEntropyFix_HartenYee(eig1, eig_lim1, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig2, eig_lim2, epsilon);
		U2NITS::Space::EigenEntropyFix_HartenYee(eig3, eig_lim3, epsilon);
	}
	// 输出特征矩阵
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
