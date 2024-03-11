#include "FluxGPU.h"
#include "../FVM_2D.h"
#include "../math/PhysicalConvertKernel.h"

void GPU::Space::Flux::calculateFluxDevice(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device) {
	// TODO: Implement this function.
    resetElementFlux(elementField_device);
    cudaDeviceSynchronize();
    calculateFlux(element_device, elementField_device, edge_device, boundary_device, infPara_device);


}

void GPU::Space::Flux::resetElementFlux(FieldSoA& elementField_device) {
    // 初始化单元数值通量
    int block_size = 512;// 最好是128 256 512
    int grid_size = (elementField_device.num + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    resetElementFluxKernel <<<grid, block>>> (elementField_device);
}

__global__ void GPU::Space::Flux::resetElementFluxKernel(FieldSoA& elementField_device) {
    // --- 核函数 --- 

    // 获取id，判断是否有效
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& i = id;
    if (i >= elementField_device.num || i < 0) return;

    // 给单元数值通量赋值
    for (int jValue = 0; jValue < 4; jValue++) {
        elementField_device.Flux[jValue][i] = 0;//单元数值通量Flux清零，为后面加减做准备
        // f->elements[ie].deltaeig = 0;//每一轮deltaeig清零 针对Roe格式 目前不需要
        // 梯度已经计算，因此不需要再次计算
    }
}

void GPU::Space::Flux::calculateFlux(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device) {

    int block_size = 512;
    int grid_size = (edge_device.num_edge + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    calculateFluxKernel <<<grid, block >>> (element_device, elementField_device, edge_device, boundary_device, infPara_device);
}

__global__ void GPU::Space::Flux::calculateFluxKernel(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device) {
    // TODO: Implement this kernel.
    // 未完成

    // 周期边界，EdgeSoA的periodicPair成员变量已完成，已经包含在edge_device中，不必再传参
    // 读取边界类型时，可以使用boundary_host

    // 获取id，判断是否有效
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& iedge = id;
    if (iedge >= edge_device.num_edge || iedge < 0) return;

    // 每条边计算无粘数值通量

    //FVM_2D* f = FVM_2D::getInstance();
    int elementL = edge_device.elementL[iedge];
    int elementR = edge_device.elementR[iedge];
    int setID = edge_device.setID[iedge];// 内部edge不属于任何set，因此setID为-1 setID的初始化见readContinueFile
    //int boundaryNum = f->boundaryManager.boundaries.size();
    int boundaryNum = boundary_device.size;
    int bType = -1;
    if (setID != -1) {
        bType = boundary_device.type[setID - 1];
    }
    double flux[4]{};

    switch (bType) {
        //对称。对欧拉方程，相当于无粘固壁
    case _BC_symmetry:
        GPU::Space::Flux::getEdgeFlux_wallNonViscous(element_device, elementField_device, edge_device, iedge, flux);
        break;
        //无粘固壁
    case _BC_wall_nonViscous:
        GPU::Space::Flux::getEdgeFlux_wallNonViscous(element_device, elementField_device, edge_device, iedge, flux);
        break;
        //无滑移绝热
    case _BC_wall_adiabat:
        // TODO: 实现无滑移绝热边界条件的计算
        GPU::Space::Flux::getEdgeFlux_wall_adiabat(element_device, elementField_device, edge_device, iedge, flux);
        break;
        //入口
    case _BC_inlet:
        GPU::Space::Flux::getEdgeFlux_farField(element_device, elementField_device, edge_device, iedge,
            infPara_device.ruvp_inlet->ptr, flux);
        break;
        //出口
    case _BC_outlet:
        GPU::Space::Flux::getEdgeFlux_farField(element_device, elementField_device, edge_device, iedge,
            infPara_device.ruvp_outlet->ptr, flux);
        break;
        //远场
    case _BC_inf:
        GPU::Space::Flux::getEdgeFlux_farField(element_device, elementField_device, edge_device, iedge,
            infPara_device.ruvp_inf->ptr, flux);
        break;
    default:// 内部：bType=-1
        if (elementR != -1) {
            // 周期和内部 统一处理
            GPU::Space::Flux::getEdgeFlux_inner_and_periodic(element_device, elementField_device, edge_device, iedge, flux);
        }
    }

    // 更新单元数值通量
    for (int j = 0; j < 4; j++) {
        elementField_device.Flux[j][elementL] += flux[j];
        // 内部边界
        if (bType == -1) {
            // 此处不是elementR==-1，防止周期边界算2次
            elementField_device.Flux[j][elementR] -= flux[j];
        }
    }
}

__device__ void GPU::Space::Flux::getEdgeFlux_wallNonViscous(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux) {
    未完待续;
}
__device__ void GPU::Space::Flux::getEdgeFlux_wall_adiabat(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux){

}
__device__ void GPU::Space::Flux::getEdgeFlux_farField(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* ruvp_inf, REAL* flux){

}
__device__ void GPU::Space::Flux::modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, const REAL* ruvp_inf){}
__device__ void GPU::Space::Flux::getEdgeFlux_inner_and_periodic(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* flux){}
__device__ void GPU::Space::Flux::getUByXYandElementID(ElementSoA& element_host, FieldSoA& elementField_host, REAL x, REAL y, int elementID, REAL* U) {}

__device__ void GPU::Space::Flux::RiemannSolve(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny, const REAL length, REAL* flux, const int scheme) {

}