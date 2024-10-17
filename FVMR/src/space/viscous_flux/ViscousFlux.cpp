#include "ViscousFlux.h"
#include "../../solvers/GPUSolver2.h"
#include "../../solvers/SolverDataGetter.h"
#include "../../global/GlobalPara.h"
#include "ViscousFluxGPU.h"

// ȡ���ҵ�Ԫ�ݶ�ƽ��ֵ
__host__ __device__ void get_U_Gradient_mean_host_device(
	myfloat& drhodx, myfloat& drhoudx, myfloat& drhovdx, myfloat& drhoEdx,
	myfloat& drhody, myfloat& drhoudy, myfloat& drhovdy, myfloat& drhoEdy,
	myint elementL, myint elementR, GPU::ElementFieldSoA& elementField_host
);

// ȡ���ҵ�Ԫ�غ���ƽ��ֵ
__host__ __device__ void get_U_mean_host_device(
	myfloat& rho, myfloat& rhou, myfloat& rhov, myfloat& rhoE,
	myint elementL, myint elementR, GPU::ElementFieldSoA& elementField_host
);

// �������������ݶȡ�������ʽ����
__host__ __device__ void get_physical_variable_and_gradient_host_device(
	myfloat& u, myfloat& v, myfloat& T,
	myfloat& dudx, myfloat& dvdx, myfloat& dTdx,
	myfloat& dudy, myfloat& dvdy, myfloat& dTdy,
	myfloat rho, myfloat rhou, myfloat rhov, myfloat rhoE,
	myfloat drhodx, myfloat drhoudx, myfloat drhovdx, myfloat drhoEdx,
	myfloat drhody, myfloat drhoudy, myfloat drhovdy, myfloat drhoEdy,
	myfloat Cv
);

// ����������ʽ(Sutherland's law)��ճ��ϵ��
myfloat get_mu_using_Sutherland_air(myfloat temperature) {
	/*
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

// ������ŵӦ��
inline void get_Reynolds_stress_2d_host_inline(
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

void viscous_flux_kernel(myint iEdge, GPU::EdgeSoA& edge_host, GPU::EdgeFieldSoA& edgeField_host, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, myfloat sutherland_C1) {
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
	get_U_Gradient_mean_host_device(
		drhodx, drhoudx, drhovdx, drhoEdx, drhody, drhoudy, drhovdy, drhoEdy,
		elementL, elementR, elementField_host
	);
	/*
	ȡ���ҵ�Ԫ�غ���ƽ��ֵ��Ϊ�����غ�����Ȼ��ת��Ϊ������
	����1����Ҫ��Roeƽ��������ƽ���Ƿ���У���ճͨ������Բ���̣���Ҫ��֤�غ�������
	����2����ʽ������ʱ�������ϴ󣬲�֪�����������Ƿ���õ�����Ҫ��������
	*/
	myfloat rho, rhou, rhov, rhoE;
	get_U_mean_host_device(rho, rhou, rhov, rhoE, elementL, elementR, elementField_host);
	myfloat u, v, T;
	myfloat dudx, dvdx, dTdx;
	myfloat dudy, dvdy, dTdy;
	get_physical_variable_and_gradient_host_device(
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
	myfloat mu = GPU::Space::get_mu_using_Sutherland_air_host_device(T, sutherland_C1);// ����ճ��ϵ��
	myfloat tau_xx_laminar, tau_yy_laminar, tau_xy_laminar;// ������ŵӦ��
	get_Reynolds_stress_2d_host_inline(
		tau_xx_laminar, tau_yy_laminar, tau_xy_laminar, mu,
		dudx, dudy, dvdx, dvdy
	);
	constexpr bool use_turbulence = false;
	myfloat mu_turbulence = 0.0;// ��ճϵ��
	myfloat tau_xx_turbulence{}, tau_yy_turbulence{}, tau_xy_turbulence{};// ������ŵӦ��
	if (use_turbulence) {
		get_Reynolds_stress_2d_host_inline(
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

	flux[0] = 0.0;
	flux[1] = (tau_xx * nx + tau_xy * ny) * area;
	flux[2] = (tau_xy * nx + tau_yy * ny) * area;
	flux[3] = (Phix * nx + Phiy * ny) * area;

	edgeField_host.Flux[0][iEdge] -= flux[0];
	edgeField_host.Flux[1][iEdge] -= flux[1];
	edgeField_host.Flux[2][iEdge] -= flux[2];
	edgeField_host.Flux[3][iEdge] -= flux[3];
	
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
	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);

	for (myint iEdge = 0; iEdge < num_edge; iEdge++) {
		viscous_flux_kernel(iEdge, edge_host, edgeField_host, element_host, elementField_host, sutherland_C1);
	}

}

