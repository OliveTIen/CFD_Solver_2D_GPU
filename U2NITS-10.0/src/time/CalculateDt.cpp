#include "CalculateDt.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"
#include "../space/viscous_flux/ViscousFlux.h"
#include "../space/viscous_flux/ViscousFluxGPU.h"
#include "../global/GlobalPara.h"


// ����edge������ͨ��alpha���Ӽ������ҵ�Ԫ��alphaC
void update_alphaC(myfloat gamma, myfloat Pr, myfloat Rcpcv, myfloat* alphaC, myfloat sutherland_C1,
	GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4], int physicsModel_equation) {

	const myfloat value_1 = (U2NITS::Math::max)(4.0 / 3.0, gamma / Pr);

	// ��ԪalphaC����
	const myint num_element = element.num_element;
	for (myint i = 0; i < num_element; i++) {
		alphaC[i] = 0.0;
	}
	// ����edge��alpha���Ӽ�����ԪalphaC
	for (myint i = 0; i < edge.num_edge; i++) {
		// ����edge�غ�������Ϊ�ڲ�edge�ͱ߽�edge�������
		myint elementL = edge.elementL[i];
		myint elementR = edge.elementR[i];
		myfloat rho, u, v, p;
		if (elementR != -1) {
			rho = 0.5 * (element_vruvp[0][elementL] + element_vruvp[0][elementR]);
			u = 0.5 * (element_vruvp[1][elementL] + element_vruvp[1][elementR]);
			v = 0.5 * (element_vruvp[2][elementL] + element_vruvp[2][elementR]);
			p = 0.5 * (element_vruvp[3][elementL] + element_vruvp[3][elementR]);
		}
		else {
			rho = element_vruvp[0][elementL];
			u = element_vruvp[1][elementL];
			v = element_vruvp[2][elementL];
			p = element_vruvp[3][elementL];
		}
		if (rho < 0 || p < 0) {
			LogWriter::logAndPrintError("Error: rho or p < 0 @U2NITS::get_global_dt_host\n");
			CExit::saveAndExit(-1);
		}

		// ����edge����ճalpha���Ӽ������൥ԪalphaC
		myfloat dx = edge.normal[0][i] * edge.length[i];
		myfloat dy = edge.normal[1][i] * edge.length[i];
		const myfloat length = edge.length[i];
		myfloat unormal = u * dx + v * dy;
		myfloat alpha = abs(unormal) + sqrt(gamma * p / rho) * sqrt(length);

		alphaC[elementL] += alpha;
		if (elementR != -1)alphaC[elementR] += alpha;

		// ����edge����ճalphaVis���Ӽ������൥ԪalphaC��ճ�Ի�����CFL�ȶ���
		if (physicsModel_equation == _EQ_NS) {
			myfloat temperature = p / Rcpcv / rho;
			myfloat mu_laminar = GPU::Space::get_mu_using_Sutherland_air_host_device(temperature, sutherland_C1);
			myfloat alphaVis = 2.0 * length * value_1 * mu_laminar / rho;

			alphaC[elementL] += alphaVis / element.volume[elementL];
			if (elementR != -1)alphaC[elementR] += alphaVis / element.volume[elementR];
		}
	}
}

// ������Ԫ��������alphaC������Ե�dt�����ȡ��Сֵ
myfloat get_minimus_dt_by_alphaC(myfloat* alphaC, const myfloat* volume, myfloat CFL, myint num_element) {
	/*
	���Ҫ�����GPU����Ļ���Ϊ��ʡ�洢�ռ䣬���Ծ���alphaC�洢dt
	�� alphaC[i] = CFL * volume[i] / alphaC[i]
	Ȼ���alphaC���й�Լ������Сֵ���Ǵ����dt
	*/
	myfloat dt = 1e10;
	for (myint i = 0; i < num_element; i++) {
		dt = (U2NITS::Math::min)(dt, CFL * volume[i] / alphaC[i]);
	}
	return dt;
}

myfloat U2NITS::Time::get_global_dt_host(myfloat currentPhysicalTime, myfloat maxPhysicalTime, myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt elementFieldVariable_dt,
	GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4], int physicsModel_equation) {
	/*
	ͳһʱ�䲽�������ڷǶ���
	�ο�CFD Pinciples and Applications P175���Լ�laminar-WangQ���� Timestep.f90

	���ȱ���edge������ͨ��alpha���Ӽ������ҵ�Ԫ��alphaC
	Ȼ�������Ԫ��������alphaC������Ե�dt�����ȡ��Сֵ
	*/

	/*
	�´��� �Ѳ��� 05-15
	*/
	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);
	update_alphaC(gamma, Pr, Rcpcv, elementFieldVariable_dt.alphaC, sutherland_C1, element, edge, element_vruvp, physicsModel_equation);
	myfloat dt = get_minimus_dt_by_alphaC(elementFieldVariable_dt.alphaC, element.volume, CFL, element.num_element);
	dt = (currentPhysicalTime + dt > maxPhysicalTime) ? (maxPhysicalTime - currentPhysicalTime) : dt;//if (t + dt > T)dt = T - t; 
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
	const myfloat value_1 = (U2NITS::Math::max)(4.0 / 3.0, gamma / Pr);
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
