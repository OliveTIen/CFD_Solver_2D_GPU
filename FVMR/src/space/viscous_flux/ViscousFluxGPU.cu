#include "ViscousFluxGPU.h"
#include "../../solvers/GPUSolver2.h"
#include "../../solvers/SolverDataGetter.h"
#include "../../global/GlobalPara.h"

// 取左右单元梯度平均值
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

	// 若elementR不存在，即是边界且非周期边界edge，则取左单元值
	if (elementR<0 || elementR>=elementField_host.num) {// has: return (iElement >= 0 && iElement < num);
		drhodx = drhodxL;
		drhoudx = drhoudxL;
		drhovdx = drhovdxL;
		drhoEdx = drhoEdxL;

		drhody = drhodyL;
		drhoudy = drhoudyL;
		drhovdy = drhovdyL;
		drhoEdy = drhoEdyL;
		return;// 用return代替else，防止分支判断，并行执行更高效
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

// 取左右单元守恒量平均值
__host__ __device__ void get_U_mean_host_device(
	myfloat& rho, myfloat& rhou, myfloat& rhov, myfloat& rhoE,
	myint elementL, myint elementR, GPU::ElementFieldSoA& elementField_host
) {
	myfloat rhoL = elementField_host.U[0][elementL];
	myfloat rhouL = elementField_host.U[1][elementL];
	myfloat rhovL = elementField_host.U[2][elementL];
	myfloat rhoEL = elementField_host.U[3][elementL];

	// 若elementR不存在，即是边界且非周期边界edge，则取左单元值
	if (elementR < 0 || elementR >= elementField_host.num) { // !elementField_host.has(elementR)
		rho = rhoL;
		rhou = rhouL;
		rhov = rhovL;
		rhoE = rhoEL;
		return;// 用return代替else
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

// 计算物理量及梯度。用求导链式法则
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
	因此
	dTdx = (dEdx - u*dudx - v*dvdx)/Cv
	*/
	myfloat one_on_Cv = 1.0 / Cv;
	myfloat V2 = u * u + v * v;
	T = one_on_Cv * (E - 0.5 * V2);
	dTdx = (dEdx - u * dudx - v * dvdx) * one_on_Cv;
	dTdy = (dEdy - u * dudy - v * dvdy) * one_on_Cv;
}

// 计算雷诺应力
__device__ inline void get_Reynolds_stress_2d_device_inline(
	myfloat& tau_xx, myfloat& tau_yy, myfloat& tau_xy,
	myfloat mu, myfloat dudx, myfloat dudy, myfloat dvdx, myfloat dvdy
) {
	// NMD: lambda/mu，第二粘性系数/粘性系数。stokes认为它等于-2/3
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
	const myfloat Cv = R / ga1;// 定容比热
	const myfloat Cp = Cv * gamma;// 定压比热

	/*
	取左右单元守恒量梯度平均值作为界面梯度。相当于梯度的常量重构，认为单元内梯度是常数
	注意这是守恒量梯度，还要转化为物理量梯度
	*/
	myfloat drhodx, drhoudx, drhovdx, drhoEdx;
	myfloat drhody, drhoudy, drhovdy, drhoEdy;
	get_U_Gradient_mean_host_device(
		drhodx, drhoudx, drhovdx, drhoEdx, drhody, drhoudy, drhovdy, drhoEdy,
		elementL, elementR, elementField_device
	);
	/*
	取左右单元守恒量平均值作为界面守恒量，然后转换为物理量
	疑问1：需要用Roe平均吗？算术平均是否可行？无粘通量是椭圆方程，需要保证守恒性质吗？
	疑问2：链式法则求导时计算量较大，不知道其他函数是否会用到？需要存起来吗？
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
	计算雷诺应力
	湍流目前先不管
	*/
	myfloat tau_xx{}, tau_yy{}, tau_xy{};// 雷诺应力。用{}初始化为0.0
	myfloat mu = GPU::Space::get_mu_using_Sutherland_air_host_device(T, sutherland_C1);// 动力粘性系数
	myfloat tau_xx_laminar, tau_yy_laminar, tau_xy_laminar;// 层流雷诺应力
	get_Reynolds_stress_2d_device_inline(
		tau_xx_laminar, tau_yy_laminar, tau_xy_laminar, mu,
		dudx, dudy, dvdx, dvdy
	);
	constexpr bool use_turbulence = false;
	myfloat mu_turbulence = 0.0;// 涡粘系数
	myfloat tau_xx_turbulence{}, tau_yy_turbulence{}, tau_xy_turbulence{};// 湍流雷诺应力
	if (use_turbulence) {
		get_Reynolds_stress_2d_device_inline(
			tau_xx_turbulence, tau_yy_turbulence, tau_xy_turbulence, mu_turbulence,
			dudx, dudy, dvdx, dvdy
		);
	}
	tau_xx = tau_xx_laminar + tau_xx_turbulence;
	tau_yy = tau_yy_laminar + tau_yy_turbulence;
	tau_xy = tau_xy_laminar + tau_xy_turbulence;

	myfloat akmu = Cp * (mu / Pr);// 传热系数。由于mu随当地温度变化，需要重新计算
	if (use_turbulence) {
		constexpr myfloat Pr_turbulence = 0.90;// 某结构程序中设置为常值0.90
		akmu += Cp * (mu_turbulence / Pr_turbulence);
	}
	/*
	计算粘性通量
	*/
	myfloat Phix = akmu * dTdx + u * tau_xx + v * tau_xy;// 粘性通量的能量项
	myfloat Phiy = akmu * dTdy + u * tau_xy + v * tau_yy;
	const myfloat nx = edge_device.normal[0][iEdge];// 边法向量(归一化)
	const myfloat ny = edge_device.normal[1][iEdge];
	const myfloat area = edge_device.length[iEdge];// 界面面积
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
	20240526 粘性通量GPU
	没有对CPU代码进行校核，直接写成GPU形式。因为加粘性后，CPU代码校核算例算得太慢了
	*/
	GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
	GPU::ElementSoA& element_device = solver->element_device;
	GPU::ElementFieldSoA& elementField_device = solver->elementField_device;
	GPU::EdgeSoA& edge_device = solver->edge_device;
	GPU::EdgeFieldSoA& edgeField_device = solver->edgeField_device;

	const myfloat R = GlobalPara::constant::R;
	constexpr myfloat Pr = 0.72;// Prantl number，某结构程序中设置为常量值0.72
	const myfloat gamma = GlobalPara::constant::gamma;
	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);

	int block_size = GPU::MY_BLOCK_SIZE;
    int grid_size = (edge_device.num_edge + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
	viscous_flux_device_kernel <<<grid, block >>> (edge_device, edgeField_device, element_device, elementField_device, R, Pr, gamma, sutherland_C1);
    getLastCudaError("edge_viscous_flux_device failed.");
}
