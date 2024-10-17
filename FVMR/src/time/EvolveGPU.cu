#include "EvolveGPU.h"
#include "../global/GlobalPara.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../space/FluxGPU.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"
#include "../space/gradient/GradientGPU.h"
#include "../math/BasicAlgorithmGPU.h"
#include <sstream>

void cuda_vector_divide_by_elements_with_weight(integer length, myfloat* dist, const myfloat* src, myfloat weight) {
    // 带权的向量对应元素相除
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
    // 计算常微分方程的右端项f=f(t,U)。与时间无关，因此简化为f(U)

    // 根据U计算Ux、Uy，即重构
    GPU::Space::Gradient::Gradient_2(element_device, node_device, edge_device, ynp_device);// 已检查
    // 根据U和Ux、Uy计算通量，存入ynp.Flux
    GPU::Space::Flux::calculateFluxDevice_2(element_device, edge_device, ynp_device);
    // 乘以体积负倒数，得到右端项f，存入ynp.Flux
    calculateFunctionF_modify_flux_by_volume(element_device.num_element, ynp_device.Flux, element_device.volume);

    getLastCudaError("calculateFunctionF_device failed.");
}

void GPU::Time::evolve_explicit_globaltimestep_device(myfloat dt, GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device) {
    myint num_element = element_device.num_element;
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_element + threadsPerBlock - 1) / threadsPerBlock;

    // 计算残值。dU/dt = f(t,U)右端项
    calculateFunctionF_device(element_device, node_device, edge_device, elementField_device);
    // 时间推进(时间积分)
    for (int i = 0; i < 4; i++) {
        GPU::Math::vector_weighted_add_kernel <<<blocksPerGrid, threadsPerBlock>>> (num_element, elementField_device.U[i], elementField_device.Flux[i], dt);
    }

    getLastCudaError("evolve_explicit_globaltimestep_device failed.");
}

void GPU::Time::evolve_explicit_localtimestep_device(GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_device) {
    myint num_element = element_device.num_element;
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_element + threadsPerBlock - 1) / threadsPerBlock;
    // 计算残值
    calculateFunctionF_device(element_device, node_device, edge_device, elementField_device);
    // 时间推进(时间积分)
    for (int i = 0; i < 4; i++) {
        GPU::Math::vector_dot_product_add_kernel <<<blocksPerGrid, threadsPerBlock>>> (num_element, elementField_device.U[i], elementField_device.Flux[i], elementFieldVariable_dt_device.alphaC);
    }

    getLastCudaError("evolve_explicit_localtimestep_device failed.");
}
