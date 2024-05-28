#include "EvolveGPU.h"
#include "../global/GlobalPara.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../space/FluxGPU.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"
#include "../space/gradient/GradientGPU.h"
#include "../math/BasicAlgorithmGPU.h"
#include <sstream>

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

    int block_size = GPU::get_max_threads_per_block();
    int grid_size = (element_device.num_element + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    TimeIntegration_1_kernel <<<grid, block >>> (dt_device, element_device, elementField_device);

}

void EvolveExplicitDevice_1(myfloat dt, GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device, GPU::EdgeSoA& edge_device, GPU::BoundarySetMap& boundary, GPU::SDevicePara& para) {
    // ��ֵͨ��
    //GPU::Space::Flux::calculateFluxDevice(element_device, elementField_device, edge_device, boundary, para);
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


void cuda_vector_divide_by_elements_with_weight(integer length, myfloat* dist, const myfloat* src, myfloat weight) {
    // ��Ȩ��������ӦԪ�����
    int threadsPerBlock = 256;
    int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;

    GPU::Math::vector_weighted_divide_kernel <<<blocksPerGrid, threadsPerBlock >>> (length, dist, src, weight);
    getLastCudaError("cuda_vector_divide_by_elements_with_weight failed.");
}

void calculateFunctionF_modify_flux_by_volume(integer num,  myfloat* ynp_flux[4], myfloat* element_volume) {

    for (int j = 0; j < 4; j++) {
        cuda_vector_divide_by_elements_with_weight(num, ynp_flux[j], element_volume, -1.0);
    }

    getLastCudaError("calculateFunctionF_modify_flux_by_volume failed.");
}

void calculateFunctionF_device(GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& ynp_device) {
    // ���㳣΢�ַ��̵��Ҷ���f=f(t,U)����ʱ���޹أ���˼�Ϊf(U)

    // ����U����Ux��Uy�����ع�
    GPU::Space::Gradient::Gradient_2(element_device, node_device, edge_device, ynp_device);// �Ѽ��
    // ����U��Ux��Uy����ͨ��������ynp.Flux
    GPU::Space::Flux::calculateFluxDevice_2(element_device, edge_device, ynp_device);
    // ����������������õ��Ҷ���f������ynp.Flux
    calculateFunctionF_modify_flux_by_volume(element_device.num_element, ynp_device.Flux, element_device.volume);

    getLastCudaError("calculateFunctionF_device failed.");
}

void cuda_vector_add_with_weight(integer length, myfloat* dist, const myfloat* src, myfloat weight) {
    // ��Ȩ�������ӷ�
    int threadsPerBlock = 256;
    int blocksPerGrid = (length + threadsPerBlock - 1) / threadsPerBlock;
    GPU::Math::vector_weighted_add_kernel <<<blocksPerGrid, threadsPerBlock >>> (length, dist, src, weight);
    getLastCudaError("cuda_vector_add_with_weight failed.");
}

void evolveSingleStep_timeIntegration_scale2Darray(integer num, myfloat* U[4], myfloat dt, myfloat* flux[4]) {
    for (int i = 0; i < 4; i++) {
        cuda_vector_add_with_weight(num, U[i], flux[i], dt);
    }

    getLastCudaError("evolveSingleStep_timeIntegration_scale2Darray failed.");
}

void GPU::Time::evolveSingleStep_device(myfloat dt, GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device) {
    // ����dU/dt = f(t,U)�Ҷ���
    calculateFunctionF_device(element_device, node_device, edge_device, elementField_device);
    // ʱ�����
    evolveSingleStep_timeIntegration_scale2Darray(elementField_device.num, elementField_device.U, dt, elementField_device.Flux);

    getLastCudaError("evolveSingleStep_device failed.");
}
