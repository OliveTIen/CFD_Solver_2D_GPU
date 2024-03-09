#include "CalculateDt.h"

// CPU�汾
REAL GPU::calculateDt(REAL t, REAL T, REAL gamma, REAL Re, REAL prl, REAL CFL, REAL Rcpcv, bool ifViscous, bool ifConstantViscous, GPU::GPUSolver2& gpuSolver) {
	// �����OpenMP��Ҫ��omp��Լ��ͣ����ܼ���ӣ���Ϊ�漰�߳̾���
	// �����GPU��ҲҪ������������
	/*
	��дCPU�汾
	���������
	prl������������Pr = Cp * mu / lambda (lambda���ȴ���ϵ��)

	*/
	REAL dt = 1e10;// a big value
	const REAL value_1 = (std::max)(4.0 / 3.0, gamma / prl);
	REAL* alpaC;
	GPU::ElementSoA& element = gpuSolver.element_host;
	GPU::EdgeSoA& edge = gpuSolver.edge_host;
	auto& element_vruvp = gpuSolver.element_vruvp;
	const long num_element = element.num_element;
	// ������Դ
	alpaC = new REAL[num_element]();// С���ű�ʾ��ʼ��Ϊ0 Ϊʲô���Ǳ���δ��ʼ����
	//for (int i = 0; i < num_element; i++) {
	//	alpaC[i] = 0.0;
	//}
	// ��ʼ����Ԫruvp ��Ԫ�غ���ת��Ԫruvp ����GPUSolver2::update_ruvp_Uold(void* _pFVM2D_);�����

	// ��ѭ�� ����laminar-WangQ Timestep.f90
	for (int i = 0; i < edge.num_edge; i++) {
		// �������غ���ֵ
		int elementL = edge.elementL[i];
		int elementR = edge.elementR[i];
		REAL rr, uu, vv, pp;
		if (elementR != -1) {
			// ��Ԫ���棬ȡ����ƽ��ֵ
			rr = 0.5 * (element_vruvp[0][elementL] + element_vruvp[0][elementR]);
			uu = 0.5 * (element_vruvp[1][elementL] + element_vruvp[1][elementR]);
			vv = 0.5 * (element_vruvp[2][elementL] + element_vruvp[2][elementR]);
			pp = 0.5 * (element_vruvp[3][elementL] + element_vruvp[3][elementR]);
		}
		else {
			// �߽��棬ȡ�ڲ൥Ԫֵ
			rr = element_vruvp[0][elementL];
			uu = element_vruvp[1][elementL];
			vv = element_vruvp[2][elementL];
			pp = element_vruvp[3][elementL];
		}
		// �������쳣ֵ������ֹ
		if (rr < 0 || pp < 0) {
			std::cout << "Error: rr or pp < 0" << std::endl;
			throw "Error: rr or pp < 0";
		}
		// ������alpha��Ȼ��ӵ����൥ԪalphaC�С���Ϊ��ճ����ճ�������
		REAL dx = edge.normal[0][i] * edge.length[i];
		REAL dy = edge.normal[1][i] * edge.length[i];
		const REAL length = edge.length[i];
		REAL unormal = uu * dx + vv * dy;
		if (ifViscous) {
			// ��ճʱ������ճ��ЧӦ������CFL�ȶ���
			// ����vmul
			REAL vmul = 0;
			if (ifConstantViscous) {
				// ��ճ��ϵ���㶨����ֱ����ճ��ϵ��
				vmul = 1.0 / Re;
			}
			else {
				REAL tt = pp / Rcpcv / rr;
				vmul = 1.458 * pow(abs(tt), 1.5) / (tt + 110.4) * 1.0 - 6;
			}
			// �������alpa+alpaVis/Volume��Ȼ��ӵ����൥Ԫ��alpaC 
			REAL alpa = abs(unormal) + sqrt(gamma * pp / rr) * sqrt(length);
			REAL alpaVis = 2.0 * (length)*value_1 * vmul / rr;
			alpaC[elementL] += alpa + alpaVis / element.volume[elementL];
			if (elementR != -1)alpaC[elementR] += alpa + alpaVis / element.volume[elementR];
		}
		else {
			// ��ճ
			REAL alpa = abs(unormal) + sqrt(gamma * pp / rr) * sqrt(length);
			alpaC[elementL] += alpa;
			if (elementR != -1)alpaC[elementR] += alpa;
		}// if ifViscous
	}// num_edge
	 // ���㵥Ԫdt
	for (int i = 0; i < num_element; i++) {
		dt = (std::min)(dt, CFL * element.volume[i] / alpaC[i]);
	}
	// ���t+dt�����趨T��������dt
	dt = (t + dt > T) ? (T - t) : dt;//if (t + dt > T)dt = T - t;
	// �ͷ���Դ
	delete[] alpaC;
	// ����dt
	return dt;
}
