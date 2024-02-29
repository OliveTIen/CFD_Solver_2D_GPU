#include "CalculateDt.h"

// CPU版本
REAL GPU::calculateDt(REAL t, REAL T, REAL gamma, REAL Re, REAL prl, REAL CFL, REAL Rcpcv, bool ifViscous, bool ifConstantViscous, GPU::GPUSolver2& gpuSolver) {
	// 如果用OpenMP，要用omp规约求和，不能简单相加，因为涉及线程竞争
	// 如果用GPU，也要考虑类似问题
	/*
	先写CPU版本
	输入参数：
	prl：普朗特数，Pr = Cp * mu / lambda (lambda：热传导系数)

	*/
	REAL dt = 1e10;// a big value
	const REAL value_1 = (std::max)(4.0 / 3.0, gamma / prl);
	REAL* alpaC;
	GPU::ElementSoA& element = gpuSolver.element_host;
	GPU::EdgeSoA& edge = gpuSolver.edge_host;
	auto& element_vruvp = gpuSolver.element_vruvp;
	const long num_element = element.num_element;
	// 申请资源
	alpaC = new REAL[num_element]();// 小括号表示初始化为0 为什么还是报错未初始化？
	//for (int i = 0; i < num_element; i++) {
	//	alpaC[i] = 0.0;
	//}
	// 初始化单元ruvp 单元守恒量转单元ruvp 已在GPUSolver2::update_ruvp_Uold(void* _pFVM2D_);中完成

	// 面循环 参照laminar-WangQ Timestep.f90
	for (int i = 0; i < edge.num_edge; i++) {
		// 计算面守恒量值
		int elementL = edge.elementL[i];
		int elementR = edge.elementR[i];
		REAL rr, uu, vv, pp;
		if (elementR != -1) {
			// 单元内面，取两侧平均值
			rr = 0.5 * (element_vruvp[0][elementL] + element_vruvp[0][elementR]);
			uu = 0.5 * (element_vruvp[1][elementL] + element_vruvp[1][elementR]);
			vv = 0.5 * (element_vruvp[2][elementL] + element_vruvp[2][elementR]);
			pp = 0.5 * (element_vruvp[3][elementL] + element_vruvp[3][elementR]);
		}
		else {
			// 边界面，取内侧单元值
			rr = element_vruvp[0][elementL];
			uu = element_vruvp[1][elementL];
			vv = element_vruvp[2][elementL];
			pp = element_vruvp[3][elementL];
		}
		// 若出现异常值，则终止
		if (rr < 0 || pp < 0) {
			std::cout << "Error: rr or pp < 0" << std::endl;
			throw "Error: rr or pp < 0";
		}
		// 计算面alpha，然后加到两侧单元alphaC中。分为有粘和无粘两种情况
		REAL dx = edge.normal[0][i] * edge.length[i];
		REAL dy = edge.normal[1][i] * edge.length[i];
		const REAL length = edge.length[i];
		REAL unormal = uu * dx + vv * dy;
		if (ifViscous) {
			// 有粘时，考虑粘性效应带来的CFL稳定性
			// 计算vmul
			REAL vmul = 0;
			if (ifConstantViscous) {
				// 若粘性系数恒定，则直接用粘性系数
				vmul = 1.0 / Re;
			}
			else {
				REAL tt = pp / Rcpcv / rr;
				vmul = 1.458 * pow(abs(tt), 1.5) / (tt + 110.4) * 1.0 - 6;
			}
			// 计算面的alpa+alpaVis/Volume，然后加到两侧单元的alpaC 
			REAL alpa = abs(unormal) + sqrt(gamma * pp / rr) * sqrt(length);
			REAL alpaVis = 2.0 * (length)*value_1 * vmul / rr;
			alpaC[elementL] += alpa + alpaVis / element.volume[elementL];
			if (elementR != -1)alpaC[elementR] += alpa + alpaVis / element.volume[elementR];
		}
		else {
			// 无粘
			REAL alpa = abs(unormal) + sqrt(gamma * pp / rr) * sqrt(length);
			alpaC[elementL] += alpa;
			if (elementR != -1)alpaC[elementR] += alpa;
		}// if ifViscous
	}// num_edge
	 // 计算单元dt
	for (int i = 0; i < num_element; i++) {
		dt = (std::min)(dt, CFL * element.volume[i] / alpaC[i]);
	}
	// 如果t+dt超出设定T，则限制dt
	dt = (t + dt > T) ? (T - t) : dt;//if (t + dt > T)dt = T - t;
	// 释放资源
	delete[] alpaC;
	// 返回dt
	return dt;
}
