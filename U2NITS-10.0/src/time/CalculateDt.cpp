#include "CalculateDt.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"


myfloat U2NITS::Time::calculateGlobalDt(myfloat currentPhysicalTime, myfloat maxPhysicalTime, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]) {
	/*
	ͳһʱ�䲽�������ڷǶ���
	�ο�CFD Pinciples and Applications P175

	*/

	bool isViscous = false;
	bool isConstantViscous = false;
	myfloat dt = 1e10;
	const myfloat value_1 = (Math::max)(4.0 / 3.0, gamma / Pr);
	const myint num_element = element.num_element;// int����ʾ21�ڣ����ȫ����myint
	// ������Դ
	myfloat* alphaC;
	alphaC = new myfloat[num_element]();
	// ��ʼ����Ԫruvp ��Ԫ�غ���ת��Ԫruvp ����GPUSolver2::update_ruvp_Uold(void* _pFVM2D_);�����

	// ��ѭ�� ����laminar-WangQ Timestep.f90
	for (myint i = 0; i < edge.num_edge; i++) {
		// �������غ���ֵ
		myint elementL = edge.elementL[i];
		myint elementR = edge.elementR[i];
		myfloat rho, u, v, p;
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
		myfloat dx = edge.normal[0][i] * edge.length[i];
		myfloat dy = edge.normal[1][i] * edge.length[i];
		const myfloat length = edge.length[i];
		myfloat unormal = u * dx + v * dy;// �����ٶȴ�С���Ա߳�
		if (isViscous) {
			// ��ճʱ������ճ��ЧӦ������CFL�ȶ���
			// ����vmul
			myfloat vmul = 0;
			if (isConstantViscous) {
				// ��ճ��ϵ���㶨����ֱ����ճ��ϵ��
				vmul = 1.0 / Re;
			}
			else {
				myfloat temperature = p / Rcpcv / rho;
				vmul = 1.458 * pow(abs(temperature), 1.5) / (temperature + 110.4) * 1.0 - 6;
			}
			// �������alpa+alpaVis/Volume��Ȼ��ӵ����൥Ԫ��alpaC 
			myfloat alpha = abs(unormal) + sqrt(gamma * p / rho) * sqrt(length);
			myfloat alphaVis = 2.0 * (length)*value_1 * vmul / rho;
			alphaC[elementL] += alpha + alphaVis / element.volume[elementL];
			if (elementR != -1)alphaC[elementR] += alpha + alphaVis / element.volume[elementR];
		}
		else {
			// ��ճ
			myfloat alpha = abs(unormal) + sqrt(gamma * p / rho) * sqrt(length);
			alphaC[elementL] += alpha;
			if (elementR != -1)alphaC[elementR] += alpha;
		}// if ifViscous
	}// num_edge
	 // ���㵥Ԫdt
	for (myint i = 0; i < num_element; i++) {
		dt = (Math::min)(dt, CFL * element.volume[i] / alphaC[i]);
	}
	// ���t+dt�����趨T��������dt
	dt = (currentPhysicalTime + dt > maxPhysicalTime) ? (maxPhysicalTime - currentPhysicalTime) : dt;//if (t + dt > T)dt = T - t;
	// �ͷ���Դ
	delete[] alphaC;
	// ����dt
	return dt;
}

void U2NITS::Time::calculateLocalTimeStep_async_Euler(myfloat& dt, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat R, myint iElement, GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]) {

	/*
	[�������ճͨ��]
	���㱾��ʱ�䲽���첽ģʽ
	�첽ģʽ����Ԫ������Ե�ͨ���������ڹ�Լ��͡��ŵ��Ǳ��ڲ��С�ռ�ÿռ�С��ȱ�������ظ�����(ÿ���汻����2��)
	��һ��ģʽ���ȼ������бߵ�ͨ������������Ȼ��ֱ�Ӽ�����Ӧ��Ԫ���ڶԵ�Ԫ�Ӽ�ʱ���漰�����ݾ���
	*/

	bool ifViscous = false;
	const myfloat _4_3 = 4.0 / 3.0;
	const myint num_element = element.num_element;

	if (!element.has(iElement)) {
		return;
	}

	myfloat lambdaConvective = 0.0;
	myfloat lambdaViscous = 0.0;
	for (int i = 0; i < 4; i++) {
		myint iEdge = element.edges[i][iElement];
		if (iEdge == -1) {
			continue;
		}
		// ����ֵ��ȡ���൥Ԫ��ֵ�����ֻ��һ�൥Ԫ����ȡһ��ֵ
		myint elementL = edge.elementL[iEdge];
		myint elementR = edge.elementR[iEdge];
		myfloat rho = element_vruvp[0][elementL];
		myfloat u = element_vruvp[1][elementL];
		myfloat v = element_vruvp[2][elementL];
		myfloat p = element_vruvp[3][elementL];
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
		myfloat nx = edge.normal[0][iEdge];
		myfloat ny = edge.normal[1][iEdge];
		myfloat normalVelocityMagnitude = Math::abs(u * nx + v * ny);
		// ����
		myfloat soundSpeed = sqrt(gamma * p / rho);
		myfloat faceArea = edge.length[iEdge];
		// ��Ԫ����ͨ���װ뾶֮��
		lambdaConvective += (normalVelocityMagnitude + soundSpeed) * faceArea;

		if (ifViscous) {// ���ա�CFDPA��P175��ʽ6.21
			throw "unimplemented";

			//myfloat a1 = Math::max(_4_3 / rho, gamma / rho);
			//myfloat a2 = 

		}
	}
	myfloat volume = element.volume[iElement];
	myfloat C = 2.0;
	dt = CFL * volume / (lambdaConvective + C * lambdaViscous);
}
