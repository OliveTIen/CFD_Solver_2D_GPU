#include "CalculateDt.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"


real U2NITS::Time::calculateGlobalDt(real currentPhysicalTime, real maxPhysicalTime, real gamma, real Re, real Pr, real CFL, real Rcpcv, GPU::ElementSoA& element, GPU::EdgeSoA& edge, real* element_vruvp[4]) {
	/*
	统一时间步长，用于非定常
	参考CFD Pinciples and Applications P175
	
	*/
	
	bool ifViscous = false;
	bool ifConstantViscous = false;
	real dt = 1e10;
	const real value_1 = (Math::max)(4.0 / 3.0, gamma / Pr);
	const integer num_element = element.num_element;// int最大表示21亿，最好全部用integer
	// 申请资源
	real* alphaC;
	alphaC = new real[num_element]();
	// 初始化单元ruvp 单元守恒量转单元ruvp 已在GPUSolver2::update_ruvp_Uold(void* _pFVM2D_);中完成

	// 面循环 参照laminar-WangQ Timestep.f90
	for (integer i = 0; i < edge.num_edge; i++) {
		// 计算面守恒量值
		integer elementL = edge.elementL[i];
		integer elementR = edge.elementR[i];
		real rho, u, v, p;
		if (elementR != -1) {
			// 单元内面，取两侧平均值
			rho = 0.5 * (element_vruvp[0][elementL] + element_vruvp[0][elementR]);
			u = 0.5 * (element_vruvp[1][elementL] + element_vruvp[1][elementR]);
			v = 0.5 * (element_vruvp[2][elementL] + element_vruvp[2][elementR]);
			p = 0.5 * (element_vruvp[3][elementL] + element_vruvp[3][elementR]);
		}
		else {
			// 边界面，取内侧单元值
			rho = element_vruvp[0][elementL];
			u = element_vruvp[1][elementL];
			v = element_vruvp[2][elementL];
			p = element_vruvp[3][elementL];
		}
		// 若出现异常值，则终止。通常是因为发散
		if (rho < 0 || p < 0) {
			LogWriter::logAndPrintError("Error: rho or p < 0 @U2NITS::calculateGlobalDt\n");
			CExit::saveAndExit(-1);
		}
		// 计算面alpha，然后加到两侧单元alphaC中。分为有粘和无粘两种情况
		real dx = edge.normal[0][i] * edge.length[i];
		real dy = edge.normal[1][i] * edge.length[i];
		const real length = edge.length[i];
		real unormal = u * dx + v * dy;// 法向速度大小乘以边长
		if (ifViscous) {
			// 有粘时，考虑粘性效应带来的CFL稳定性
			// 计算vmul
			real vmul = 0;
			if (ifConstantViscous) {
				// 若粘性系数恒定，则直接用粘性系数
				vmul = 1.0 / Re;
			}
			else {
				real temperature = p / Rcpcv / rho;
				vmul = 1.458 * pow(abs(temperature), 1.5) / (temperature + 110.4) * 1.0 - 6;
			}
			// 计算面的alpa+alpaVis/Volume，然后加到两侧单元的alpaC 
			real alpha = abs(unormal) + sqrt(gamma * p / rho) * sqrt(length);
			real alphaVis = 2.0 * (length)*value_1 * vmul / rho;
			alphaC[elementL] += alpha + alphaVis / element.volume[elementL];
			if (elementR != -1)alphaC[elementR] += alpha + alphaVis / element.volume[elementR];
		}
		else {
			// 无粘
			real alpha = abs(unormal) + sqrt(gamma * p / rho) * sqrt(length);
			alphaC[elementL] += alpha;
			if (elementR != -1)alphaC[elementR] += alpha;
		}// if ifViscous
	}// num_edge
	 // 计算单元dt
	for (integer i = 0; i < num_element; i++) {
		dt = (Math::min)(dt, CFL * element.volume[i] / alphaC[i]);
	}
	// 如果t+dt超出设定T，则限制dt
	dt = (currentPhysicalTime + dt > maxPhysicalTime) ? (maxPhysicalTime - currentPhysicalTime) : dt;//if (t + dt > T)dt = T - t;
	// 释放资源
	delete[] alphaC;
	// 返回dt
	return dt;
}

void U2NITS::Time::calculateLocalTimeStep_async_Euler(real& dt, real gamma, real Re, real Pr, real CFL, real R, integer iElement, GPU::ElementSoA& element, GPU::EdgeSoA& edge, real* element_vruvp[4]) {
	
	/*
	[仅完成无粘通量]
	计算本地时间步，异步模式
	异步模式：单元计算各自的通量，不存在规约求和。优点是便于并行、占用空间小，缺点是有重复计算(每个面被计算2遍)
	另一种模式是先计算所有边的通量，存起来，然后分别加减到对应单元。在对单元加减时，涉及到数据竞争
	*/

	bool ifViscous = false;
	const real _4_3 = 4.0 / 3.0;
	const integer num_element = element.num_element;

	if (!element.has(iElement)) {
		return;
	}

	real lambdaConvective = 0.0;
	real lambdaViscous = 0.0;
	for (int i = 0; i < 4; i++) {
		integer iEdge = element.edges[i][iElement];
		if (iEdge == -1) {
			continue;
		}
		// 界面值，取两侧单元均值。如果只有一侧单元，则取一侧值
		integer elementL = edge.elementL[iEdge];
		integer elementR = edge.elementR[iEdge];
		real rho = element_vruvp[0][elementL];
		real u = element_vruvp[1][elementL];
		real v = element_vruvp[2][elementL];
		real p = element_vruvp[3][elementL];
		if (element.has(elementR)) {
			rho += element_vruvp[0][elementR];
			u += element_vruvp[1][elementR];
			v += element_vruvp[2][elementR];
			p += element_vruvp[3][elementR];
			rho *= 0.5;
			u *= 0.5;
			v *= 0.5;
			p *= 0.5;
		}
		// 法向速度大小
		real nx = edge.normal[0][iEdge];
		real ny = edge.normal[1][iEdge];
		real normalVelocityMagnitude = Math::abs(u * nx + v * ny);
		// 声速
		real soundSpeed = sqrt(gamma * p / rho);
		real faceArea = edge.length[iEdge];
		// 单元对流通量谱半径之和
		lambdaConvective += (normalVelocityMagnitude + soundSpeed) * faceArea;

		if (ifViscous) {// 参照《CFDPA》P175公式6.21
			throw "unimplemented";

			//real a1 = Math::max(_4_3 / rho, gamma / rho);
			//real a2 = 

		}
	}
	real volume = element.volume[iElement];
	real C = 2.0;
	dt = CFL * volume / (lambdaConvective + C * lambdaViscous);
}
