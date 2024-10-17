#include "ViscousFluxGPU.h"
#include "../../solvers/GPUSolver2.h"
#include "../../solvers/SolverDataGetter.h"
#include "../../global/GlobalPara.h"

// ȡ���ҵ�Ԫ�ݶ�ƽ��ֵ
__host__ __device__ void get_U_Gradient_mean_host_device(
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
	if (elementR<0 || elementR>=elementField_host.num) {// has: return (iElement >= 0 && iElement < num);
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
__host__ __device__ void get_U_mean_host_device(
	myfloat& rho, myfloat& rhou, myfloat& rhov, myfloat& rhoE,
	myint elementL, myint elementR, GPU::ElementFieldSoA& elementField_host
) {
	myfloat rhoL = elementField_host.U[0][elementL];
	myfloat rhouL = elementField_host.U[1][elementL];
	myfloat rhovL = elementField_host.U[2][elementL];
	myfloat rhoEL = elementField_host.U[3][elementL];

	// ��elementR�����ڣ����Ǳ߽��ҷ����ڱ߽�edge����ȡ��Ԫֵ
	if (elementR < 0 || elementR >= elementField_host.num) { // !elementField_host.has(elementR)
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
__host__ __device__ void get_physical_variable_and_gradient_host_device(
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

// ������ŵӦ��
__device__ inline void get_Reynolds_stress_2d_device_inline(
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

__global__ void viscous_flux_device_kernel(GPU::EdgeSoA edge_device, GPU::EdgeFieldSoA edgeField_device, GPU::ElementSoA element_device, GPU::ElementFieldSoA elementField_device
	,myfloat R, myfloat Pr, myfloat gamma, myfloat sutherland_C1
) {
	const myint iEdge = blockIdx.x * blockDim.x + threadIdx.x;
	if (iEdge >= edge_device.num_edge || iEdge < 0) return;

	const myint elementL = edge_device.elementL[iEdge];
	const myint elementR = edge_device.elementR[iEdge];
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
		elementL, elementR, elementField_device
	);
	/*
	ȡ���ҵ�Ԫ�غ���ƽ��ֵ��Ϊ�����غ�����Ȼ��ת��Ϊ������
	����1����Ҫ��Roeƽ��������ƽ���Ƿ���У���ճͨ������Բ���̣���Ҫ��֤�غ�������
	����2����ʽ������ʱ�������ϴ󣬲�֪�����������Ƿ���õ�����Ҫ��������
	*/
	myfloat rho, rhou, rhov, rhoE;
	get_U_mean_host_device(rho, rhou, rhov, rhoE, elementL, elementR, elementField_device);
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
	get_Reynolds_stress_2d_device_inline(
		tau_xx_laminar, tau_yy_laminar, tau_xy_laminar, mu,
		dudx, dudy, dvdx, dvdy
	);
	constexpr bool use_turbulence = false;
	myfloat mu_turbulence = 0.0;// ��ճϵ��
	myfloat tau_xx_turbulence{}, tau_yy_turbulence{}, tau_xy_turbulence{};// ������ŵӦ��
	if (use_turbulence) {
		get_Reynolds_stress_2d_device_inline(
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
		akmu += Cp * (mu_turbulence / Pr_turbulence);
	}
	/*
	����ճ��ͨ��
	*/
	myfloat Phix = akmu * dTdx + u * tau_xx + v * tau_xy;// ճ��ͨ����������
	myfloat Phiy = akmu * dTdy + u * tau_xy + v * tau_yy;
	const myfloat nx = edge_device.normal[0][iEdge];// �߷�����(��һ��)
	const myfloat ny = edge_device.normal[1][iEdge];
	const myfloat area = edge_device.length[iEdge];// �������
	myfloat flux[4]{};

	flux[0] = 0.0;
	flux[1] = (tau_xx * nx + tau_xy * ny) * area;
	flux[2] = (tau_xy * nx + tau_yy * ny) * area;
	flux[3] = (Phix * nx + Phiy * ny) * area;

	edgeField_device.Flux[0][iEdge] -= flux[0];
	edgeField_device.Flux[1][iEdge] -= flux[1];
	edgeField_device.Flux[2][iEdge] -= flux[2];
	edgeField_device.Flux[3][iEdge] -= flux[3];
	
}

void GPU::Space::edge_viscous_flux_device() {
	/*
	20240526 ճ��ͨ��GPU
	û�ж�CPU�������У�ˣ�ֱ��д��GPU��ʽ����Ϊ��ճ�Ժ�CPU����У���������̫����
	*/
	GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
	GPU::ElementSoA& element_device = solver->element_device;
	GPU::ElementFieldSoA& elementField_device = solver->elementField_device;
	GPU::EdgeSoA& edge_device = solver->edge_device;
	GPU::EdgeFieldSoA& edgeField_device = solver->edgeField_device;

	const myfloat R = GlobalPara::constant::R;
	constexpr myfloat Pr = 0.72;// Prantl number��UNITs������Ϊ����ֵ0.72
	const myfloat gamma = GlobalPara::constant::gamma;
	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);

	int block_size = GPU::MY_BLOCK_SIZE;
    int grid_size = (edge_device.num_edge + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
	viscous_flux_device_kernel <<<grid, block >>> (edge_device, edgeField_device, element_device, elementField_device, R, Pr, gamma, sutherland_C1);
    getLastCudaError("edge_viscous_flux_device failed.");
}
