#include "EvolveGPU.h"
#include "../global/GlobalPara.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../space/FluxGPU.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"

void GPU::Time::EvolveDevice(REAL dt, int flag_timeAdvance, ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary, DGlobalPara& para) {
    if (flag_timeAdvance == _EVO_explicit) {
        EvolveExplicitDevice(dt, element_device, elementField_device, edge_device, boundary, para);
    }
    else {
        LogWriter::logAndPrintError("Error: invalid evolve method.\n");
        exit(1);
    }
}

void GPU::Time::EvolveExplicitDevice(REAL dt, ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary, DGlobalPara& para) {
    // 数值通量
    GPU::Space::Flux::calculateFluxDevice(element_device, elementField_device, edge_device, boundary, para);
    cudaDeviceSynchronize();
    catchCudaErrorAndExit();

    // 时间积分
    TimeIntegration(dt, element_device, elementField_device);

}

void GPU::Time::TimeIntegration(REAL dt, ElementSoA& element_device, FieldSoA& elementField_device) {
    DReal dt_device(&dt);

    int block_size = 512;// 最好是128 256 512
    int grid_size = (element_device.num_element + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    TimeIntegrationKernel <<<grid, block >>> (dt_device, element_device, elementField_device);

}

__global__ void GPU::Time::TimeIntegrationKernel(DReal& dt, ElementSoA& element, FieldSoA& elementField) {
    // 获取id，判断是否有效
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& ie = id;
    if (ie >= element.num_element || ie < 0) return;

    REAL omega = element.volume[ie];
    for (int j = 0; j < 4; j++) {
        elementField.U[j][ie] -= *(dt.ptr) / omega * elementField.Flux[j][ie];
    }

}
