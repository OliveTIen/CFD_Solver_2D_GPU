#include "ViscousFlux.h"
#include "../../solvers/GPUSolver2.h"
#include "../../solvers/SolverDataGetter.h"
#include "../../global/GlobalPara.h"
#include "ViscousFluxGPU.h"

// 取左右单元梯度平均值
__host__ __device__ void get_U_Gradient_mean_host_device(
	myfloat& drhodx, myfloat& drhoudx, myfloat& drhovdx, myfloat& drhoEdx,
	myfloat& drhody, myfloat& drhoudy, myfloat& drhovdy, myfloat& drhoEdy,
	myint elementL, myint elementR, GPU::ElementFieldSoA& elementField_host
);

// 取左右单元守恒量平均值
__host__ __device__ void get_U_mean_host_device(
	myfloat& rho, myfloat& rhou, myfloat& rhov, myfloat& rhoE,
	myint elementL, myint elementR, GPU::ElementFieldSoA& elementField_host
);

// 计算物理量及梯度。用求导链式法则
__host__ __device__ void get_physical_variable_and_gradient_host_device(
	myfloat& u, myfloat& v, myfloat& T,
	myfloat& dudx, myfloat& dvdx, myfloat& dTdx,
	myfloat& dudy, myfloat& dvdy, myfloat& dTdy,
	myfloat rho, myfloat rhou, myfloat rhov, myfloat rhoE,
	myfloat drhodx, myfloat drhoudx, myfloat drhovdx, myfloat drhoEdx,
	myfloat drhody, myfloat drhoudy, myfloat drhovdy, myfloat drhoEdy,
	myfloat Cv
);

// 用苏世兰公式(Sutherland's law)求粘性系数
myfloat get_mu_using_Sutherland_air(myfloat temperature) {
	/*
	用静态变量，第一次计算，后续无需重复计算
	静态变量生存期为程序周期。下次调用时仍使用上次的值
	问题：
	某些情况下inf rho=0导致发散
	*/
	static myfloat C_sutherland = 0.0;
	static bool is_called_first_time = true;
	if (is_called_first_time) {
		myfloat T_inf = 690.0;// 远场边界条件的无穷远静温 某original.org中取690
		C_sutherland = 117.0 / T_inf;
		is_called_first_time = false;
	}
	myfloat& T = temperature;
	return (1.0 + C_sutherland) / (T + C_sutherland) * pow(T, 1.5);
}

// 计算雷诺应力
inline void get_Reynolds_stress_2d_host_inline(
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

void viscous_flux_kernel(myint iEdge, GPU::EdgeSoA& edge_host, GPU::EdgeFieldSoA& edgeField_host, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, myfloat sutherland_C1) {
	/*
	对iEdge，计算数值通量，并
	创建：tgl, 20240430
	校核：
	*/
	const myint elementL = edge_host.elementL[iEdge];
	const myint elementR = edge_host.elementR[iEdge];
	const myfloat R = GlobalPara::constant::R;
	constexpr myfloat Pr = 0.72;// Prantl number，某结构程序中设置为常量值0.72
	const myfloat gamma = GlobalPara::constant::gamma;
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
		elementL, elementR, elementField_host
	);
	/*
	取左右单元守恒量平均值作为界面守恒量，然后转换为物理量
	疑问1：需要用Roe平均吗？算术平均是否可行？无粘通量是椭圆方程，需要保证守恒性质吗？
	疑问2：链式法则求导时计算量较大，不知道其他函数是否会用到？需要存起来吗？
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
	计算雷诺应力
	湍流目前先不管
	*/
	myfloat tau_xx{}, tau_yy{}, tau_xy{};// 雷诺应力。用{}初始化为0.0
	myfloat mu = GPU::Space::get_mu_using_Sutherland_air_host_device(T, sutherland_C1);// 动力粘性系数
	myfloat tau_xx_laminar, tau_yy_laminar, tau_xy_laminar;// 层流雷诺应力
	get_Reynolds_stress_2d_host_inline(
		tau_xx_laminar, tau_yy_laminar, tau_xy_laminar, mu,
		dudx, dudy, dvdx, dvdy
	);
	constexpr bool use_turbulence = false;
	myfloat mu_turbulence = 0.0;// 涡粘系数
	myfloat tau_xx_turbulence{}, tau_yy_turbulence{}, tau_xy_turbulence{};// 湍流雷诺应力
	if (use_turbulence) {
		get_Reynolds_stress_2d_host_inline(
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
		akmu += Cp * (mu_turbulence/ Pr_turbulence);
	}
	/*
	计算粘性通量。此处并没有用到雷诺数
	
	雷诺数是对方程无量纲化后的产物，一般不需要写进方程
	某结构程序依照远场参数进行了归一化，因此用了雷诺数

	无粘通量量纲(->表示量纲相同)
	F[1]   -> rho*u*u+p          -> [Pa]
	F[3]   -> (rho*E+p)*u        -> [Pa * m/s]

	粘性通量量纲
	F_v[1] -> tau_xx -> 2*mu*u   -> [Pa]
	F_v[3] -> Phix   -> u*tau_xx -> [Pa * m/s]

	可以看到粘性通量与无粘通量量纲相同。对于数值通量，还需乘以面积，因为要用到Gauss公式
	*/
	myfloat Phix = akmu * dTdx + u * tau_xx + v * tau_xy;// 粘性通量的能量项
	myfloat Phiy = akmu * dTdy + u * tau_xy + v * tau_yy;
	const myfloat nx = edge_host.normal[0][iEdge];// 边法向量(归一化)
	const myfloat ny = edge_host.normal[1][iEdge];
	const myfloat area = edge_host.length[iEdge];// 界面面积
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
	20240430 粘性通量
	计算粘性通量后，对edgeField_host的值进行修改
	将循环体提取成函数viscous_flux_kernel，与CUDA kernel保持一致，便于转化
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

