#include "ViscousFlux.h"
#include "../../solvers/GPUSolver2.h"
#include "../../solvers/SolverDataGetter.h"
#include "../../global/GlobalPara.h"

// ȡ���ҵ�Ԫ�ݶ�ƽ��ֵ
void get_U_Gradient_mean(
	myfloat& drhodx, myfloat& drhoudx, myfloat& drhovdx, myfloat& drhoEdx,
	myfloat& drhody, myfloat& drhoudy, myfloat& drhovdy, myfloat& drhoEdy,
	myint elementL, myint elementR, GPU::ElementFieldSoA& elementField_host
) {
	myfloat drhodxL = elementField_host.Ux[0][elementL];
	myfloat drhoudxL = elementField_host.Ux[1][elementL];
	myfloat drhovdxL = elementField_host.Ux[2][elementL];
	myfloat drhoEdxL = elementField_host.Ux[3][elementL];

	myfloat drhodyL = elementField_host.Uy[0][elementL];
	myfloat drhoudyL = elementField_host.Uy[1][elementL];
	myfloat drhovdyL = elementField_host.Uy[2][elementL];
	myfloat drhoEdyL = elementField_host.Uy[3][elementL];

	// ��elementR�����ڣ����Ǳ߽��ҷ����ڱ߽�edge����ȡ��Ԫֵ
	if (!elementField_host.has(elementR)) {
		drhodx = drhodxL;
		drhoudx = drhoudxL;
		drhovdx = drhovdxL;
		drhoEdx = drhoEdxL;

		drhody = drhodyL;
		drhoudy = drhoudyL;
		drhovdy = drhovdyL;
		drhoEdy = drhoEdyL;
		return;// ��return����else����ֹ��֧�жϣ�����ִ�и���Ч
	}

	myfloat drhodxR = elementField_host.Ux[0][elementR];
	myfloat drhoudxR = elementField_host.Ux[1][elementR];
	myfloat drhovdxR = elementField_host.Ux[2][elementR];
	myfloat drhoEdxR = elementField_host.Ux[3][elementR];

	myfloat drhodyR = elementField_host.Uy[0][elementR];
	myfloat drhoudyR = elementField_host.Uy[1][elementR];
	myfloat drhovdyR = elementField_host.Uy[2][elementR];
	myfloat drhoEdyR = elementField_host.Uy[3][elementR];

	drhodx = 0.5 * (drhodxL + drhodxR);
	drhoudx = 0.5 * (drhoudxL + drhoudxR);
	drhovdx = 0.5 * (drhovdxL + drhovdxR);
	drhoEdx = 0.5 * (drhoEdxL + drhoEdxR);

	drhody = 0.5 * (drhodyL + drhodyR);
	drhoudy = 0.5 * (drhoudyL + drhoudyR);
	drhovdy = 0.5 * (drhovdyL + drhovdyR);
	drhoEdy = 0.5 * (drhoEdyL + drhoEdyR);

}

// ȡ���ҵ�Ԫ�غ���ƽ��ֵ
void get_U_mean(
	myfloat& rho, myfloat& rhou, myfloat& rhov, myfloat& rhoE,
	myint elementL, myint elementR, GPU::ElementFieldSoA& elementField_host
) {
	myfloat rhoL = elementField_host.U[0][elementL];
	myfloat rhouL = elementField_host.U[1][elementL];
	myfloat rhovL = elementField_host.U[2][elementL];
	myfloat rhoEL = elementField_host.U[3][elementL];

	// ��elementR�����ڣ����Ǳ߽��ҷ����ڱ߽�edge����ȡ��Ԫֵ
	if (!elementField_host.has(elementR)) {
		rho = rhoL;
		rhou = rhouL;
		rhov = rhovL;
		rhoE = rhoEL;
		return;// ��return����else
	}

	myfloat rhoR = elementField_host.U[0][elementR];
	myfloat rhouR = elementField_host.U[1][elementR];
	myfloat rhovR = elementField_host.U[2][elementR];
	myfloat rhoER = elementField_host.U[3][elementR];

	rho = 0.5 * (rhoL + rhoR);
	rhou = 0.5 * (rhouL + rhouR);
	rhov = 0.5 * (rhovL + rhovR);
	rhoE = 0.5 * (rhoEL + rhoER);
	
}

// �������������ݶȡ�������ʽ����
void get_physical_variable_and_gradient(
	myfloat& u, myfloat& v, myfloat& T,
	myfloat& dudx, myfloat& dvdx, myfloat& dTdx,
	myfloat& dudy, myfloat& dvdy, myfloat& dTdy,
	myfloat rho, myfloat rhou, myfloat rhov, myfloat rhoE,
	myfloat drhodx, myfloat drhoudx, myfloat drhovdx, myfloat drhoEdx,
	myfloat drhody, myfloat drhoudy, myfloat drhovdy, myfloat drhoEdy,
	myfloat Cv
) {
	myfloat one_on_rho = 1.0 / rho;
	u = rhou * one_on_rho;
	v = rhov * one_on_rho;
	myfloat E = rhoE * one_on_rho;


	dudx = one_on_rho * (drhoudx - u * drhodx);
	dvdx = one_on_rho * (drhovdx - v * drhodx);
	dudy = one_on_rho * (drhoudy - u * drhody);
	dvdy = one_on_rho * (drhovdy - v * drhody);

	myfloat dEdx = one_on_rho * (drhoEdx - E * drhodx);
	myfloat dEdy = one_on_rho * (drhoEdy - E * drhody);

	/*
	E = e + 0.5 * V2 = Cv * T + 0.5 * (u*u+v*v)
	dEdx = Cv * dTdx + 0.5 * d(u*u+v*v)dx = Cv * dTdx + u*dudx + v*dvdx
	���
	dTdx = (dEdx - u*dudx - v*dvdx)/Cv
	*/
	myfloat one_on_Cv = 1.0 / Cv;
	myfloat V2 = u * u + v * v;
	T = one_on_Cv * (E - 0.5 * V2);
	dTdx = (dEdx - u * dudx - v * dvdx) * one_on_Cv;
	dTdy = (dEdy - u * dudy - v * dvdy) * one_on_Cv;
}

// ����������ʽ(Sutherland's law)��ճ��ϵ��
myfloat get_mu_using_Sutherland_air(myfloat temperature) {
	/*

	����UNITs, Csthlnd = 117.0/T_inf
	�þ�̬��������һ�μ��㣬���������ظ�����
	��̬����������Ϊ�������ڡ��´ε���ʱ��ʹ���ϴε�ֵ
	���⣺
	ĳЩ�����inf rho=0���·�ɢ
	*/
	static myfloat C_sutherland = 0.0;
	static bool is_called_first_time = true;
	if (is_called_first_time) {
		myfloat T_inf = 690.0;// Զ���߽�����������Զ���£�UNITs original.org��ȡ690
		C_sutherland = 117.0 / T_inf;
		is_called_first_time = false;
	}
	myfloat& T = temperature;
	return (1.0 + C_sutherland) / (T + C_sutherland) * pow(T, 1.5);
}

myfloat get_mu_using_Sutherland_air_2(myfloat temperature) {
	/*
	�˴��������ڿ�������

	��������������� https://www.cfd-online.com/Wiki/Sutherland%27s_law
		mu0 = 1.716e-5;// �ο�
		T0 = 273.15;// �ο��¶�
		S = 110.4;
		C1 = mu0 / pow(T0, 1.5) * (T0 + S) = 1.45793e-06������վ����Ǻ�
	*/
	myfloat& T = temperature;
	myfloat S = 110.4;
	myfloat C1 = 1.45793e-06;
	return C1 * pow(T, 1.5) / (T + S);
}

// ������ŵӦ��
inline void get_Reynolds_stress_2d(
	myfloat& tau_xx, myfloat& tau_yy, myfloat& tau_xy,
	myfloat mu, myfloat dudx, myfloat dudy, myfloat dvdx, myfloat dvdy
) {
	// NMD: lambda/mu���ڶ�ճ��ϵ��/ճ��ϵ����stokes��Ϊ������-2/3
	constexpr myfloat NMD = -2.0 / 3.0;
	myfloat div = dudx + dvdy;
	tau_xx = mu * (NMD * div + 2.0 * dudx);
	tau_yy = mu * (NMD * div + 2.0 * dvdy);
	tau_xy = mu * (dudy + dvdx);
}

void viscous_flux_kernel(myint iEdge, GPU::EdgeSoA& edge_host, GPU::EdgeFieldSoA& edgeField_host, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host) {
	/*
	��iEdge��������ֵͨ������
	������tgl, 20240430
	У�ˣ�
	*/
	const myint elementL = edge_host.elementL[iEdge];
	const myint elementR = edge_host.elementR[iEdge];
	const myfloat R = GlobalPara::constant::R;
	constexpr myfloat Pr = 0.72;// Prantl number��UNITs������Ϊ����ֵ0.72
	const myfloat gamma = GlobalPara::constant::gamma;
	const myfloat ga1 = gamma - 1.0;
	const myfloat Cv = R / ga1;// ���ݱ���
	const myfloat Cp = Cv * gamma;// ��ѹ����
	/*
	ȡ���ҵ�Ԫ�غ����ݶ�ƽ��ֵ��Ϊ�����ݶȡ��൱���ݶȵĳ����ع�����Ϊ��Ԫ���ݶ��ǳ���
	ע�������غ����ݶȣ���Ҫת��Ϊ�������ݶ�
	*/
	myfloat drhodx, drhoudx, drhovdx, drhoEdx;
	myfloat drhody, drhoudy, drhovdy, drhoEdy;
	get_U_Gradient_mean(
		drhodx, drhoudx, drhovdx, drhoEdx, drhody, drhoudy, drhovdy, drhoEdy,
		elementL, elementR, elementField_host
	);
	/*
	ȡ���ҵ�Ԫ�غ���ƽ��ֵ��Ϊ�����غ�����Ȼ��ת��Ϊ������
	����1����Ҫ��Roeƽ��������ƽ���Ƿ���У���ճͨ������Բ���̣���Ҫ��֤�غ�������
	����2����ʽ������ʱ�������ϴ󣬲�֪�����������Ƿ���õ�����Ҫ��������
	*/
	myfloat rho, rhou, rhov, rhoE;
	get_U_mean(rho, rhou, rhov, rhoE, elementL, elementR, elementField_host);
	myfloat u, v, T;
	myfloat dudx, dvdx, dTdx;
	myfloat dudy, dvdy, dTdy;
	get_physical_variable_and_gradient(
		u, v, T,
		dudx, dvdx, dTdx,
		dudy, dvdy, dTdy,
		rho, rhou, rhov, rhoE,
		drhodx, drhoudx, drhovdx, drhoEdx,
		drhody, drhoudy, drhovdy, drhoEdy,
		Cv
	);
	/*
	������ŵӦ��
	����Ŀǰ�Ȳ���
	*/
	myfloat tau_xx{}, tau_yy{}, tau_xy{};// ��ŵӦ������{}��ʼ��Ϊ0.0
	myfloat mu = get_mu_using_Sutherland_air_2(T);// ����ճ��ϵ��
	myfloat tau_xx_laminar, tau_yy_laminar, tau_xy_laminar;// ������ŵӦ��
	get_Reynolds_stress_2d(
		tau_xx_laminar, tau_yy_laminar, tau_xy_laminar, mu,
		dudx, dudy, dvdx, dvdy
	);
	constexpr bool use_turbulence = false;
	myfloat mu_turbulence = 0.0;// ��ճϵ��
	myfloat tau_xx_turbulence{}, tau_yy_turbulence{}, tau_xy_turbulence{};// ������ŵӦ��
	if (use_turbulence) {
		get_Reynolds_stress_2d(
			tau_xx_turbulence, tau_yy_turbulence, tau_xy_turbulence, mu_turbulence,
			dudx, dudy, dvdx, dvdy
		);
	}
	tau_xx = tau_xx_laminar + tau_xx_turbulence;
	tau_yy = tau_yy_laminar + tau_yy_turbulence;
	tau_xy = tau_xy_laminar + tau_xy_turbulence;

	myfloat akmu = Cp * (mu / Pr);// ����ϵ��������mu�浱���¶ȱ仯����Ҫ���¼���
	if (use_turbulence) {
		constexpr myfloat Pr_turbulence = 0.90;// UNITs������Ϊ��ֵ0.90
		akmu += Cp * (mu_turbulence/ Pr_turbulence);
	}
	/*
	����ճ��ͨ�����˴���û���õ���ŵ��
	
	��ŵ���ǶԷ��������ٻ���Ĳ��һ�㲻��Ҫд������
	UNITs����Զ�����������˹�һ�������������ŵ��

	��ճͨ������(->��ʾ������ͬ)
	F[1]   -> rho*u*u+p          -> [Pa]
	F[3]   -> (rho*E+p)*u        -> [Pa * m/s]

	ճ��ͨ������
	F_v[1] -> tau_xx -> 2*mu*u   -> [Pa]
	F_v[3] -> Phix   -> u*tau_xx -> [Pa * m/s]

	���Կ���ճ��ͨ������ճͨ��������ͬ��������ֵͨ������������������ΪҪ�õ�Gauss��ʽ
	*/
	myfloat Phix = akmu * dTdx + u * tau_xx + v * tau_xy;// ճ��ͨ����������
	myfloat Phiy = akmu * dTdy + u * tau_xy + v * tau_yy;
	const myfloat nx = edge_host.normal[0][iEdge];// �߷�����(��һ��)
	const myfloat ny = edge_host.normal[1][iEdge];
	const myfloat area = edge_host.length[iEdge];// �������
	myfloat flux[4]{};
	//myfloat Re = GlobalPara::constant::Re;// ��һ���������е�λ1����Ӧ�������ȳ߶ȵ���ŵ��
	//myfloat area_on_Re = area / Re;
	flux[0] = 0.0;
	flux[1] = (tau_xx * nx + tau_xy * ny) * area;
	flux[2] = (tau_xy * nx + tau_yy * ny) * area;
	flux[3] = (Phix * nx + Phiy * ny) * area;


	edgeField_host.Flux[0][iEdge] -= flux[0];
	edgeField_host.Flux[1][iEdge] -= flux[1];
	edgeField_host.Flux[2][iEdge] -= flux[2];
	edgeField_host.Flux[3][iEdge] -= flux[3];
	// �����Ǽӻ��Ǽ������ҵ�Ԫ���������ҪУ�ˣ��������Ƿ���ȷ
}




void U2NITS::Space::edge_viscous_flux() {
	/*
	20240430 ճ��ͨ��
	����ճ��ͨ���󣬶�edgeField_host��ֵ�����޸�

	��ѭ������ȡ�ɺ���viscous_flux_kernel����CUDA kernel����һ�£�����ת��
	*/

	GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
	GPU::ElementSoA& element_host = solver->element_host;
	GPU::ElementFieldSoA& elementField_host = solver->elementField_host;
	GPU::EdgeSoA& edge_host = solver->edge_host;
	GPU::EdgeFieldSoA& edgeField_host = solver->edgeField_host;

	myint num_edge = edge_host.num_edge;

	for (myint iEdge = 0; iEdge < num_edge; iEdge++) {
		viscous_flux_kernel(iEdge, edge_host, edgeField_host, element_host, elementField_host);
	}

}

