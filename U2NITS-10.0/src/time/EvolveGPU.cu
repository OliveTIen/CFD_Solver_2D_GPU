#include "EvolveGPU.h"
#include "../global/GlobalPara.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../space/FluxGPU.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"
#include "../space/gradient/GradientGPU.h"
#include <sstream>

void checkCudaErrorAndExit(const char* file, int line, cudaError_t err) {
    if (err != cudaSuccess) {
        std::stringstream ss;
        ss << "failed at " << file << ", line " << line << ", error code " << cudaGetErrorString(err) << "\n";
        LogWriter::logAndPrintError(ss.str());
        exit(EXIT_FAILURE);
    }
}

__global__ void TimeIntegration_1_kernel(GPU::DReal& dt, GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField) {
    // ��ȡid���ж��Ƿ���Ч
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& ie = id;
    if (ie >= element.num_element || ie < 0) return;

    myfloat omega = element.volume[ie];
    for (int j = 0; j < 4; j++) {
        elementField.U[j][ie] -= *(dt.ptr) / omega * elementField.Flux[j][ie];
    }

}

void TimeIntegration_1(myfloat dt, GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device) {
    GPU::DReal dt_device(&dt);

    int block_size = 512;// �����128 256 512
    int grid_size = (element_device.num_element + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    TimeIntegration_1_kernel << <grid, block >> > (dt_device, element_device, elementField_device);

}

void EvolveExplicitDevice_1(myfloat dt, GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device, GPU::EdgeSoA& edge_device, GPU::BoundarySetMap& boundary, GPU::SDevicePara& para) {
    // ��ֵͨ��
    GPU::Space::Flux::calculateFluxDevice(element_device, elementField_device, edge_device, boundary, para);
    cudaDeviceSynchronize();
    GPU::catchCudaErrorAndExit();

    // ʱ�����
    TimeIntegration_1(dt, element_device, elementField_device);

}

void GPU::Time::EvolveDevice(myfloat dt, int flag_timeAdvance, ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary, SDevicePara& para) {
    if (flag_timeAdvance == _EVO_explicit) {
        EvolveExplicitDevice_1(dt, element_device, elementField_device, edge_device, boundary, para);
    }
    else {
        LogWriter::logAndPrintError("Error: invalid evolve method.\n");
        exit(1);
    }
}

__global__ void cuda_vector_divide_by_elements_with_weight_kernel(integer length, myfloat* dist, const myfloat* src, myfloat weight) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < length) {
        dist[i] /= weight * src[i];// �������ȼ� �����Ҳ�˷�������*=��/=
    }
}

void cuda_vector_divide_by_elements_with_weight(integer length, myfloat* dist, const myfloat* src, myfloat weight) {
    // ��Ȩ��������ӦԪ�����
    int threadsPerBlock = 256;
    int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
    cuda_vector_divide_by_elements_with_weight_kernel <<<blocksPerGrid, threadsPerBlock >>> (length, dist, src, weight);
    checkCudaErrorAndExit(__FILE__, __LINE__, cudaGetLastError());
}

void calculateFunctionF_modify_flux_by_volume(integer num,  myfloat* ynp_flux[4], myfloat* element_volume) {

    for (int j = 0; j < 4; j++) {
        //for (integer ie = 0; ie < num; ie++) {
        //    myfloat minus_one_on_volume = -1.0 / element_volume[ie];
        //    ynp_flux[j][ie] = minus_one_on_volume * ynp_flux[j][ie];
        //}
        cuda_vector_divide_by_elements_with_weight(num, ynp_flux[j], element_volume, -1.0);
    }
}

void calculateFunctionF(GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& ynp, GPU::BoundarySetMap& boundary, GPU::SDevicePara& para) {
    // ���㳣΢�ַ��̵��Ҷ���f=f(t,U)����ʱ���޹أ���˼�Ϊf(U)

    // ����U����Ux��Uy�����ع�
    GPU::Space::Gradient::Gradient_2(element_device, ynp, edge_device);
    // ����U��Ux��Uy����ͨ��������ynp.Flux
    GPU::Space::Flux::calculateFluxDevice_2(element_device, ynp, edge_device, boundary, para);
    // ����������������õ��Ҷ���f������ynp.Flux
    calculateFunctionF_modify_flux_by_volume(element_device.num_element, ynp.Flux, element_device.volume);

/*
  dim3 grid((dx / TILEX) + (!(dx % TILEX) ? 0 : 1),
            (dy / TILEY) + (!(dy % TILEY) ? 0 : 1));
  dim3 tids(TIDSX, TIDSY);

  updateVelocity_k<<<grid, tids>>>(v, vx, vy, dx, pdx, dy, TILEY / TIDSY,
                                   tPitch);
  getLastCudaError("updateVelocity_k failed.");
*/
}

__global__ void cuda_vector_add_with_weight_kernel(integer length, myfloat* dist, const myfloat* src, myfloat weight) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < length) {
        dist[i] += weight * src[i];
    }
}

void cuda_vector_add_with_weight(integer length, myfloat* dist, const myfloat* src, myfloat weight) {
    // ��Ȩ�������ӷ�
    int threadsPerBlock = 256;
    int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
    cuda_vector_add_with_weight_kernel <<<blocksPerGrid, threadsPerBlock >>> (length, dist, src, weight);
    checkCudaErrorAndExit(__FILE__, __LINE__, cudaGetLastError());
}

void evolveSingleStep_timeIntegration_scale2Darray(integer num, myfloat* U[4], myfloat dt, myfloat* flux[4]) {
    for (int i = 0; i < 4; i++) {
        cuda_vector_add_with_weight(num, U[i], flux[i], dt);
    }
}

void GPU::Time::evolveSingleStep_device(myfloat dt, GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device, BoundarySetMap& boundary, SDevicePara& para) {

    // ����dU/dt = f(t,U)�Ҷ���
    calculateFunctionF(element_device, node_device, edge_device, elementField_device, boundary, para);
    // ʱ�����
    evolveSingleStep_timeIntegration_scale2Darray(elementField_device.num, elementField_device.U, dt, elementField_device.Flux);
}
