#include "CalculateDt.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"


real U2NITS::Time::calculateGlobalDt(real currentPhysicalTime, real maxPhysicalTime, real gamma, real Re, real Pr, real CFL, real Rcpcv, GPU::ElementSoA& element, GPU::EdgeSoA& edge, real* element_vruvp[4]) {
	/*
	ͳһʱ�䲽�������ڷǶ���
	�ο�CFD Pinciples and Applications P175
	
	*/
	
	bool ifViscous = false;
	bool ifConstantViscous = false;
	real dt = 1e10;
	const real value_1 = (Math::max)(4.0 / 3.0, gamma / Pr);
	const integer num_element = element.num_element;// int����ʾ21�ڣ����ȫ����integer
	// ������Դ
	real* alphaC;
	alphaC = new real[num_element]();
	// ��ʼ����Ԫruvp ��Ԫ�غ���ת��Ԫruvp ����GPUSolver2::update_ruvp_Uold(void* _pFVM2D_);�����

	// ��ѭ�� ����laminar-WangQ Timestep.f90
	for (integer i = 0; i < edge.num_edge; i++) {
		// �������غ���ֵ
		integer elementL = edge.elementL[i];
		integer elementR = edge.elementR[i];
		real rho, u, v, p;
		if (elementR != -1) {
			// ��Ԫ���棬ȡ����ƽ��ֵ
			rho = 0.5 * (element_vruvp[0][elementL] + element_vruvp[0][elementR]);
			u = 0.5 * (element_vruvp[1][elementL] + element_vruvp[1][elementR]);
			v = 0.5 * (element_vruvp[2][elementL] + element_vruvp[2][elementR]);
			p = 0.5 * (element_vruvp[3][elementL] + element_vruvp[3][elementR]);
		}
		else {
			// �߽��棬ȡ�ڲ൥Ԫֵ
			rho = element_vruvp[0][elementL];
			u = element_vruvp[1][elementL];
			v = element_vruvp[2][elementL];
			p = element_vruvp[3][elementL];
		}
		// �������쳣ֵ������ֹ��ͨ������Ϊ��ɢ
		if (rho < 0 || p < 0) {
			LogWriter::logAndPrintError("Error: rho or p < 0 @U2NITS::calculateGlobalDt\n");
			CExit::saveAndExit(-1);
		}
		// ������alpha��Ȼ��ӵ����൥ԪalphaC�С���Ϊ��ճ����ճ�������
		real dx = edge.normal[0][i] * edge.length[i];
		real dy = edge.normal[1][i] * edge.length[i];
		const real length = edge.length[i];
		real unormal = u * dx + v * dy;// �����ٶȴ�С���Ա߳�
		if (ifViscous) {
			// ��ճʱ������ճ��ЧӦ������CFL�ȶ���
			// ����vmul
			real vmul = 0;
			if (ifConstantViscous) {
				// ��ճ��ϵ���㶨����ֱ����ճ��ϵ��
				vmul = 1.0 / Re;
			}
			else {
				real temperature = p / Rcpcv / rho;
				vmul = 1.458 * pow(abs(temperature), 1.5) / (temperature + 110.4) * 1.0 - 6;
			}
			// �������alpa+alpaVis/Volume��Ȼ��ӵ����൥Ԫ��alpaC 
			real alpha = abs(unormal) + sqrt(gamma * p / rho) * sqrt(length);
			real alphaVis = 2.0 * (length)*value_1 * vmul / rho;
			alphaC[elementL] += alpha + alphaVis / element.volume[elementL];
			if (elementR != -1)alphaC[elementR] += alpha + alphaVis / element.volume[elementR];
		}
		else {
			// ��ճ
			real alpha = abs(unormal) + sqrt(gamma * p / rho) * sqrt(length);
			alphaC[elementL] += alpha;
			if (elementR != -1)alphaC[elementR] += alpha;
		}// if ifViscous
	}// num_edge
	 // ���㵥Ԫdt
	for (integer i = 0; i < num_element; i++) {
		dt = (Math::min)(dt, CFL * element.volume[i] / alphaC[i]);
	}
	// ���t+dt�����趨T��������dt
	dt = (currentPhysicalTime + dt > maxPhysicalTime) ? (maxPhysicalTime - currentPhysicalTime) : dt;//if (t + dt > T)dt = T - t;
	// �ͷ���Դ
	delete[] alphaC;
	// ����dt
	return dt;
}

void U2NITS::Time::calculateLocalTimeStep_async_Euler(real& dt, real gamma, real Re, real Pr, real CFL, real R, integer iElement, GPU::ElementSoA& element, GPU::EdgeSoA& edge, real* element_vruvp[4]) {
	
	/*
	[�������ճͨ��]
	���㱾��ʱ�䲽���첽ģʽ
	�첽ģʽ����Ԫ������Ե�ͨ���������ڹ�Լ��͡��ŵ��Ǳ��ڲ��С�ռ�ÿռ�С��ȱ�������ظ�����(ÿ���汻����2��)
	��һ��ģʽ���ȼ������бߵ�ͨ������������Ȼ��ֱ�Ӽ�����Ӧ��Ԫ���ڶԵ�Ԫ�Ӽ�ʱ���漰�����ݾ���
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
		// ����ֵ��ȡ���൥Ԫ��ֵ�����ֻ��һ�൥Ԫ����ȡһ��ֵ
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
		// �����ٶȴ�С
		real nx = edge.normal[0][iEdge];
		real ny = edge.normal[1][iEdge];
		real normalVelocityMagnitude = Math::abs(u * nx + v * ny);
		// ����
		real soundSpeed = sqrt(gamma * p / rho);
		real faceArea = edge.length[iEdge];
		// ��Ԫ����ͨ���װ뾶֮��
		lambdaConvective += (normalVelocityMagnitude + soundSpeed) * faceArea;

		if (ifViscous) {// ���ա�CFDPA��P175��ʽ6.21
			throw "unimplemented";

			//real a1 = Math::max(_4_3 / rho, gamma / rho);
			//real a2 = 

		}
	}
	real volume = element.volume[iElement];
	real C = 2.0;
	dt = CFL * volume / (lambdaConvective + C * lambdaViscous);
}
