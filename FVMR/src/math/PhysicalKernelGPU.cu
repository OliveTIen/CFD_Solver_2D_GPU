#include "PhysicalKernelGPU.h"
#include "../gpu/GPUGlobalFunction.h"


__global__ void gpu_math_get_ruvp_by_U_device_kernel(myfloat* ruvp[4], const myfloat* const U[4], myint length, myfloat gamma) {
    // U:rho,rho_u,rho_v,rho_E ruvp:rho,u,v,p
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < length) {
        ruvp[0][i] = U[0][i];
        ruvp[1][i] = U[1][i] / U[0][i];
        ruvp[2][i] = U[2][i] / U[0][i];
        ruvp[3][i] = ruvp[0][i] * (gamma - 1) * (U[3][i] / U[0][i] - (ruvp[1][i] * ruvp[1][i] + ruvp[2][i] * ruvp[2][i]) * 0.5);
    }

}

void get_ruvp_by_U_device(myfloat* ruvp_device[4], const myfloat* const U_device[4], myint length, myfloat gamma) {
	int block_size = GPU::get_max_threads_per_block();
	//int block_size = GPU::MY_BLOCK_SIZE;
	int grid_size = (length + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	gpu_math_get_ruvp_by_U_device_kernel <<<grid, block>>> (ruvp_device, U_device, length, gamma);
	getLastCudaError("get_ruvp_by_U_device failed.");
}

__global__ void gpu_math_get_ruvp_by_U_device_kernel_2(GPU::ElementFieldSoA elementField_device, myfloat gamma) {
    // U:rho,rho_u,rho_v,rho_E ruvp:rho,u,v,p
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    auto& ruvp = elementField_device.ruvp;
    auto& U = elementField_device.U;
    if (i < elementField_device.num) {
        ruvp[0][i] = U[0][i];
        ruvp[1][i] = U[1][i] / U[0][i];
        ruvp[2][i] = U[2][i] / U[0][i];
        ruvp[3][i] = ruvp[0][i] * (gamma - 1) * (U[3][i] / U[0][i] - (ruvp[1][i] * ruvp[1][i] + ruvp[2][i] * ruvp[2][i]) * 0.5);
    }

}

void get_ruvp_by_U_device_2(GPU::ElementFieldSoA& elementField_device, myfloat gamma) {
    int block_size = GPU::get_max_threads_per_block();
    //int block_size = GPU::MY_BLOCK_SIZE;
    int grid_size = (elementField_device.num + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    gpu_math_get_ruvp_by_U_device_kernel_2 <<<grid, block>>>(elementField_device, gamma);
    getLastCudaError("get_ruvp_by_U_device failed.");
}


void GPU::Math::update_ruvp_Uold_device_2(GPU::ElementFieldSoA& elementField_device, myfloat gamma) {
    //cudaDeviceSynchronize();
    for (int i = 0; i < 4; i++) {
        cudaMemcpy(elementField_device.Uold[i], elementField_device.U[i], elementField_device.num * sizeof(myfloat), cudaMemcpyDeviceToDevice);
    }
    getLastCudaError("cudaMemcpy device to device failed.");

    get_ruvp_by_U_device_2(elementField_device, gamma);
}
