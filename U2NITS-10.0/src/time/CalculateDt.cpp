#include "CalculateDt.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"
#include "../space/viscous_flux/ViscousFlux.h"
#include "../space/viscous_flux/ViscousFluxGPU.h"
#include "../global/GlobalPara.h"


// ��ԪalphaC���㣬Ȼ�����edge��alpha���Ӽ�����ԪalphaC
void update_alphaC_host(myfloat gamma, myfloat Pr, myfloat Rcpcv, myfloat* alphaC, myfloat sutherland_C1,
	GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation) {

	const myfloat value_1 = (U2NITS::Math::max)(4.0 / 3.0, gamma / Pr);
	auto& element_vruvp = elementField.ruvp;

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

// �õ�Ԫ���Ե�alphaC����dt���������Դ���alphaC����
void get_local_dt_using_alphaC_CFL_volume_host(myfloat* alphaC, const myfloat* volume, myfloat CFL, myint array_size) {
	for (myint i = 0; i < array_size; i++) {
		alphaC[i] = CFL * volume[i] / alphaC[i];
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

myfloat U2NITS::Time::get_global_dt_host(myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_host,
	GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation) {
	/*
	ͳһʱ�䲽�������ڷǶ���
	�ο�CFD Pinciples and Applications P175���Լ�laminar-WangQ���� Timestep.f90

	���ȱ���edge������ͨ��alpha���Ӽ������ҵ�Ԫ��alphaC
	Ȼ�������Ԫ��������alphaC������Ե�dt�����ȡ��Сֵ
	*/

	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);
	update_alphaC_host(gamma, Pr, Rcpcv, elementFieldVariable_dt_host.alphaC, sutherland_C1, element, edge, elementField, physicsModel_equation);
	myfloat dt = get_minimus_dt_by_alphaC(elementFieldVariable_dt_host.alphaC, element.volume, CFL, element.num_element);

	return dt;
}

myfloat U2NITS::Time::get_global_dt_host_0529beta(myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_host, GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation) {
	// ������ beta�棬��GPU����һ��
	
	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);
	update_alphaC_host(gamma, Pr, Rcpcv, elementFieldVariable_dt_host.alphaC, sutherland_C1, element, edge, elementField, physicsModel_equation);
	get_local_dt_using_alphaC_CFL_volume_host(elementFieldVariable_dt_host.alphaC, element.volume, CFL, element.num_element);
	myfloat dt = U2NITS::Math::get_min_of_vector(elementFieldVariable_dt_host.alphaC, element.num_element);

	return dt;
}

void U2NITS::Time::calculateLocalTimeStep_async_Euler_kernel(myfloat& dt, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat R, myint iElement, GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]) {

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

void U2NITS::Time::get_local_dt_host(myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_host, GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation) {
	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);
	update_alphaC_host(gamma, Pr, Rcpcv, elementFieldVariable_dt_host.alphaC, sutherland_C1, element, edge, elementField, physicsModel_equation);
	get_local_dt_using_alphaC_CFL_volume_host(elementFieldVariable_dt_host.alphaC, element.volume, CFL, element.num_element);
}
