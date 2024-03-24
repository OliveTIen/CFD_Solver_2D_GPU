#include "FluxGPU.h"
#include "../FVM_2D.h"
#include "../math/Math.h"
#include "convection/ConvectionGPU.h"

void GPU::Space::Flux::calculateFluxDevice(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DGlobalPara& infPara_device) {
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

void GPU::Space::Flux::calculateFlux(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DGlobalPara& infPara_device) {

    int block_size = 512;
    int grid_size = (edge_device.num_edge + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    calculateFluxKernel <<<grid, block >>> (element_device, elementField_device, edge_device, boundary_device, infPara_device);
}

__global__ void GPU::Space::Flux::calculateFluxKernel(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, BoundarySetMap& boundary_device, DGlobalPara& para) {
    // TODO: Implement this kernel.
    // 未完成

    // 周期边界，EdgeSoA的periodicPair成员变量已完成，已经包含在edge_device中，不必再传参
    // 读取边界类型时，可以使用boundary_host

    // 获取id，判断是否有效
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& iedge = id;
    if (iedge >= edge.num_edge || iedge < 0) return;

    // 每条边计算无粘数值通量

    //FVM_2D* f = FVM_2D::getInstance();
    int elementL = edge.elementL[iedge];
    int elementR = edge.elementR[iedge];
    int setID = edge.setID[iedge];// 内部edge不属于任何set，因此setID为-1 setID的初始化见readContinueFile
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
        GPU::Space::Flux::getEdgeFlux_wallNonViscous(element, elementField, edge, iedge, flux, para);
        break;
        //无粘固壁
    case _BC_wall_nonViscous:
        GPU::Space::Flux::getEdgeFlux_wallNonViscous(element, elementField, edge, iedge, flux, para);
        break;
        //无滑移绝热
    case _BC_wall_adiabat:
        // TODO: 实现无滑移绝热边界条件的计算
        GPU::Space::Flux::getEdgeFlux_wall_adiabat(element, elementField, edge, iedge, flux, para);
        break;
        //入口
    case _BC_inlet:
        GPU::Space::Flux::getEdgeFlux_farField(element, elementField, edge, iedge, flux, para);
        break;
        //出口
    case _BC_outlet:
        GPU::Space::Flux::getEdgeFlux_farField(element, elementField, edge, iedge, flux, para);
        break;
        //远场
    case _BC_inf:
        GPU::Space::Flux::getEdgeFlux_farField(element, elementField, edge, iedge, flux, para);
        break;
    default:// 内部：bType=-1
        if (elementR != -1) {
            // 周期和内部 统一处理
            GPU::Space::Flux::getEdgeFlux_inner_and_periodic(element, elementField, edge, iedge, flux, para);
        }
    }

    // 更新单元数值通量
    for (int j = 0; j < 4; j++) {
        elementField.Flux[j][elementL] += flux[j];
        // 内部边界
        if (bType == -1) {
            // 此处不是elementR==-1，防止周期边界算2次
            elementField.Flux[j][elementR] -= flux[j];
        }
    }
}

__device__ void GPU::Space::Flux::getEdgeFlux_wallNonViscous(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, long idx, REAL* flux, DGlobalPara& para) {
    //未完待续;
    // 函数中idx相当于iedge
    // 该函数分别计算边的左右两侧U的极限，然后求解flux

    REAL nx = edge.normal[0][idx];// normal已初始化
    REAL ny = edge.normal[1][idx];
    REAL x_edge = edge.xy[0][idx];
    REAL y_edge = edge.xy[1][idx];
    REAL length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // 计算U_L，即edge U的左极限
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementL, U_L);

    // 计算U_R 用对称性
    REAL U_R[4]{};
    REAL uL = U_L[1] / U_L[0];
    REAL vL = U_L[2] / U_L[0];
    REAL uR = uL - 2 * (uL * nx + vL * ny) * nx;
    REAL vR = vL - 2 * (uL * nx + vL * ny) * ny;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    // 计算flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::getEdgeFlux_wall_adiabat(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, long idx, REAL* flux, DGlobalPara& infPara_device){

}
__device__ void GPU::Space::Flux::getEdgeFlux_farField(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, long idx, REAL* flux, DGlobalPara& infPara_device){

}
__device__ void GPU::Space::Flux::modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, DGlobalPara& infPara_device){}
__device__ void GPU::Space::Flux::getEdgeFlux_inner_and_periodic(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, long idx, REAL* flux, DGlobalPara& infPara_device){}
__device__ void GPU::Space::Flux::getUByXYandElementID(ElementSoA& element_device, FieldSoA& elementField_device, REAL x, REAL y, int elementID, REAL* U) {}

__device__ void GPU::Space::Flux::RiemannSolve(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny,
    const REAL length, REAL* flux, DGlobalPara& para) {
    real faceNormal[2]{ nx,ny };

    switch (*(para.scheme->ptr)) {
    case _SOL_LocalLaxFriedrichs:
        GPU::Space::Convection::LocalLaxFriedrichs2d(UL, UR, nx, ny, length, flux, *(para.constant_gamma->ptr));
        break;
    case _SOL_Roe:
        GPU::Space::Convection::ConvectRoeCommon2d(UL, UR, faceNormal, length, flux, para);
        break;
    default:
        break;
    }
}