#include "CalculateDt.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"
#include "../space/viscous_flux/ViscousFlux.h"
#include "../space/viscous_flux/ViscousFluxGPU.h"
#include "../global/GlobalPara.h"
#include "../math/Math.h"

// 单元alphaC归零，然后计算edge的alpha，加减到单元alphaC
void update_alphaC_host(myfloat gamma, myfloat Pr, myfloat Rcpcv, myfloat* alphaC, myfloat sutherland_C1,
	GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation) {

	const myfloat value_1 = (U2NITS::Math::max)(4.0 / 3.0, gamma / Pr);
	auto& element_vruvp = elementField.ruvp;

	// 单元alphaC归零
	const myint num_element = element.num_element;
	for (myint i = 0; i < num_element; i++) {
		alphaC[i] = 0.0;
	}
	// 计算edge的alpha，加减到单元alphaC
	for (myint i = 0; i < edge.num_edge; i++) {
		// 计算edge守恒量，分为内部edge和边界edge两种情况
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

		// 计算edge的无粘alpha，加减到两侧单元alphaC
		myfloat dx = edge.normal[0][i] * edge.length[i];
		myfloat dy = edge.normal[1][i] * edge.length[i];
		const myfloat length = edge.length[i];
		myfloat unormal = u * dx + v * dy;
		myfloat alpha = abs(unormal) + sqrt(gamma * p / rho) * sqrt(length);

		alphaC[elementL] += alpha;
		if (elementR != -1)alphaC[elementR] += alpha;

		// 计算edge的有粘alphaVis，加减到两侧单元alphaC。粘性会增加CFL稳定性
		if (physicsModel_equation == _EQ_NS) {
			myfloat temperature = p / Rcpcv / rho;
			myfloat mu_laminar = GPU::Space::get_mu_using_Sutherland_air_host_device(temperature, sutherland_C1);
			myfloat alphaVis = 2.0 * length * value_1 * mu_laminar / rho;

			alphaC[elementL] += alphaVis / element.volume[elementL];
			if (elementR != -1)alphaC[elementR] += alphaVis / element.volume[elementR];
		}
	}
}

// 用单元各自的alphaC计算dt，计算结果仍存入alphaC数组
void get_local_dt_using_alphaC_CFL_volume_host(myfloat* alphaC, const myfloat* volume, myfloat CFL, myint array_size) {
	for (myint i = 0; i < array_size; i++) {
		alphaC[i] = CFL * volume[i] / alphaC[i];
	}
}

// 遍历单元，各自用alphaC计算各自的dt，最后取最小值
myfloat get_minimus_dt_by_alphaC(myfloat* alphaC, const myfloat* volume, myfloat CFL, myint num_element) {
	/*
	如果要改造成GPU代码的话，为节省存储空间，可以就用alphaC存储dt
	即 alphaC[i] = CFL * volume[i] / alphaC[i]
	然后对alphaC进行规约，其最小值就是待求的dt
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
	统一时间步长，用于非定常
	参考CFD Pinciples and Applications P175，以及laminar-WangQ代码 Timestep.f90

	首先遍历edge，计算通量alpha，加减到左右单元的alphaC
	然后遍历单元，各自用alphaC计算各自的dt，最后取最小值
	*/

	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);
	update_alphaC_host(gamma, Pr, Rcpcv, elementFieldVariable_dt_host.alphaC, sutherland_C1, element, edge, elementField, physicsModel_equation);
	myfloat dt = get_minimus_dt_by_alphaC(elementFieldVariable_dt_host.alphaC, element.volume, CFL, element.num_element);

	return dt;
}

myfloat U2NITS::Time::get_global_dt_host_0529beta(myfloat gamma, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_host, GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation) {
	// 测试中 beta版，跟GPU保持一致
	
	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);
	update_alphaC_host(gamma, Pr, Rcpcv, elementFieldVariable_dt_host.alphaC, sutherland_C1, element, edge, elementField, physicsModel_equation);
	get_local_dt_using_alphaC_CFL_volume_host(elementFieldVariable_dt_host.alphaC, element.volume, CFL, element.num_element);
	myfloat dt = U2NITS::Math::get_min_of_vector(elementFieldVariable_dt_host.alphaC, element.num_element);

	return dt;
}

void U2NITS::Time::calculateLocalTimeStep_async_Euler_kernel(myfloat& dt, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat R, myint iElement, GPU::ElementSoA& element, GPU::EdgeSoA& edge, myfloat* element_vruvp[4]) {

	/*
	[仅完成无粘通量]
	计算本地时间步，异步模式
	异步模式：单元计算各自的通量，不存在规约求和。优点是便于并行、占用空间小，缺点是有重复计算(每个面被计算2遍)
	另一种模式是先计算所有边的通量，存起来，然后分别加减到对应单元。在对单元加减时，涉及到数据竞争
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
		// 界面值，取两侧单元均值。如果只有一侧单元，则取一侧值
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
		// 法向速度大小
		myfloat nx = edge.normal[0][iEdge];
		myfloat ny = edge.normal[1][iEdge];
		myfloat normalVelocityMagnitude = Math::abs(u * nx + v * ny);
		// 声速
		myfloat soundSpeed = sqrt(gamma * p / rho);
		myfloat faceArea = edge.length[iEdge];
		// 单元对流通量谱半径之和
		lambdaConvective += (normalVelocityMagnitude + soundSpeed) * faceArea;

		if (ifViscous) {// 参照《CFDPA》P175公式6.21
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
