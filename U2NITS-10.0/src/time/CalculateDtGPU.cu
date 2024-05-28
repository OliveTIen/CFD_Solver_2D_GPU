#include "CalculateDtGPU.h"
#include "../output/LogWriter.h"
#include "../math/MathGPU.h"
#include "../space/viscous_flux/ViscousFluxGPU.h"
#include "../global/GlobalPara.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../math/ReduceGPU.h"

// update_alphaC_device的核函数
__global__ void update_alphaC_device_kernel(myfloat gamma, myfloat Pr, myfloat Rcpcv, myfloat* alphaC, myfloat sutherland_C1,
	GPU::ElementSoA element, GPU::EdgeSoA edge, GPU::ElementFieldSoA elementField, int physicsModel_equation) {

	// 获取id，判断是否有效
	const myint iEdge = blockIdx.x * blockDim.x + threadIdx.x;
	if (iEdge >= edge.num_edge || iEdge < 0) return;
	auto& i = iEdge;
	auto& element_vruvp = elementField.ruvp;

	// 计算edge守恒量，分为内部edge和边界edge两种情况
	const myfloat value_1 = (GPU::Math::max)(4.0 / 3.0, gamma / Pr);
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
		printf("Error: rho or p < 0 @U2NITS::update_alphaC_device_kernel\n");
		//exit(-1); do nothing
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

// 单元alphaC归零，然后计算edge的alpha，加减到单元alphaC
void update_alphaC_device(myfloat gamma, myfloat Pr, myfloat Rcpcv, myfloat* alphaC, myfloat sutherland_C1, 
	GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation) {

	const myfloat value_1 = (GPU::Math::max)(4.0 / 3.0, gamma / Pr);

	// 单元alphaC归零
	cudaMemset(alphaC, 0, element.num_element * sizeof(myfloat));
	getLastCudaError("clear_element_alphaC_device failed.");

	// 计算edge的alpha，加减到单元alphaC
	int block_size = GPU::get_max_threads_per_block();
	int grid_size = (edge.num_edge + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	update_alphaC_device_kernel <<<grid, block>>> (gamma, Pr, Rcpcv, alphaC, sutherland_C1, element, edge, elementField, physicsModel_equation);
	getLastCudaError("update_alphaC_device failed.");
	
}

// 用单元各自的alphaC计算dt，并存入alphaC数组
void get_local_dt_using_alphaC_CFL_volume(myfloat* alphaC, const myfloat* volume, myfloat CFL, myint length) {
	/*
	计算 dt ，把结果存入alphaC数组，即alphaC = CFL*volume/alphaC，使用kernel: “带权的向量倒数”
	*/
	int block_size = GPU::get_max_threads_per_block();
	int grid_size = (length + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	GPU::Math::vector_weighted_reciprocal_kernel <<<grid, block>>> (length, alphaC, volume, CFL);
	getLastCudaError("get_local_dt_using_alphaC_CFL_volume failed.");
}

void GPU::Time::get_global_dt_device(myfloat currentPhysicalTime, myfloat maxPhysicalTime, myfloat gamma, myfloat Re, myfloat Pr, myfloat CFL, myfloat Rcpcv, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_device, GPU::ElementSoA& element, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, int physicsModel_equation) {
	/*
	涉及到规约操作。注意reduce_device中n必须为2的幂，用num_reduce而不是num_element
	*/

	//auto& ccc = elementFieldVariable_dt_device.alphaC;

	auto& element_var = elementFieldVariable_dt_device;
	myfloat sutherland_C1 = GPU::Space::get_Sutherland_C1_host_device(GPU::Space::S_Sutherland, GlobalPara::constant::mu0, GlobalPara::constant::T0);
	update_alphaC_device(gamma, Pr, Rcpcv, element_var.alphaC, sutherland_C1, element, edge, elementField, physicsModel_equation);
	get_local_dt_using_alphaC_CFL_volume(element_var.alphaC, element.volume, CFL, element.num_element);
	GPU::Math::reduce_device(element_var.num_reduce, element_var.alphaC, element_var.dev_output, false, GPU::Math::ReduceType::reduceType_min);
	getLastCudaError("get_global_dt_device failed.");
}