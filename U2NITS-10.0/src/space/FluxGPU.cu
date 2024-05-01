#include "FluxGPU.h"
#include "../FVM_2D.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../math/MathGPU.h"
#include "convection/ConvectionGPU.h"
#include "../output/LogWriter.h"
#include "restrict/RestrictGPU.h"

__global__ void resetElementFluxKernel_1(GPU::ElementFieldSoA elementField_device) {
    // --- 核函数 --- 

    // 获取id，判断是否有效
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& i = id;
    if (i >= elementField_device.num || i < 0) return;

    // 给单元数值通量赋值
    for (int jValue = 0; jValue < 4; jValue++) {
        elementField_device.Flux[jValue][i] = 0;//单元数值通量Flux清零，为后面加减做准备 Invalid __global__ write of size 8 bytes
        // f->elements[ie].deltaeig = 0;//每一轮deltaeig清零 针对Roe格式 目前不需要
        // 梯度已经计算，因此不需要再次计算
    }
}

void resetElementFlux_1(GPU::ElementFieldSoA& elementField_device) {
    // 初始化单元数值通量
    int block_size = 512;// 最好是128 256 512
    int grid_size = (elementField_device.num + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    resetElementFluxKernel_1 <<<grid, block>>> (elementField_device);
}

__device__ void GPU::Space::Flux::getEdgeFlux_wallNonViscous_kernel(ElementSoA& element, ElementFieldSoA& elementField, EdgeSoA& edge, long idx, myfloat* flux, SDevicePara para) {
    //未完待续;
    // 函数中idx相当于iedge
    // 该函数分别计算边的左右两侧U的极限，然后求解flux

    myfloat nx = edge.normal[0][idx];// normal已初始化
    myfloat ny = edge.normal[1][idx];
    myfloat x_edge = edge.xy[0][idx];
    myfloat y_edge = edge.xy[1][idx];
    myfloat length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // 计算U_L，即edge U的左极限
    myfloat U_L[4]{};
    Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementL, U_L, para);

    // 计算U_R 用对称性
    myfloat U_R[4]{};
    myfloat uL = U_L[1] / U_L[0];
    myfloat vL = U_L[2] / U_L[0];
    myfloat uR = uL - 2 * (uL * nx + vL * ny) * nx;
    myfloat vR = vL - 2 * (uL * nx + vL * ny) * ny;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    // 计算flux
    RiemannSolve_kernel(U_L, U_R, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::getEdgeFlux_wall_adiabat_kernel(ElementSoA& element, ElementFieldSoA& elementField, EdgeSoA& edge, long idx, myfloat* flux, SDevicePara para) {
    // 固壁边界。函数中idx相当于iedge

    myfloat nx = edge.normal[0][idx];// normal已初始化
    myfloat ny = edge.normal[1][idx];
    myfloat x_edge = edge.xy[0][idx];
    myfloat y_edge = edge.xy[1][idx];
    myfloat length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // 计算U_L，即edge U的左极限
    myfloat U_L[4]{};
    Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementL, U_L, para);

    // 计算U_R 用对称性
    myfloat U_R[4]{};
    myfloat uL = U_L[1] / U_L[0];
    myfloat vL = U_L[2] / U_L[0];
    // 法向坐标系下 UnL = uL * nx + vL * ny，VnL = -uL * ny + vL * nx
    // UnR = - UnL，VnR = -VnL
    // 原坐标系下 uR = unR*nx-vnR*ny， vR=unR*ny+vnR*nx
    // 因此 uR = -uL 还不如直接取反

    // 按UNITS中BC_WALL.f90的算法，速度取反，可是结果仍然不变
    // 不对 结果云图中没有体现粘性
    //myfloat uR = uL - 2 * (uL * nx + vL * ny) * nx;
    //myfloat vR = vL - 2 * (uL * nx + vL * ny) * ny;
    myfloat uR = -uL;
    myfloat vR = -vL;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    // 计算flux
    RiemannSolve_kernel(U_L, U_R, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::getEdgeFlux_farField_kernel(ElementSoA& element, ElementFieldSoA& elementField, EdgeSoA& edge, long idx, myfloat* flux, SDevicePara para) {
    // 函数功能：根据远场边界条件计算边界数值通量 已完成

    myfloat nx = edge.normal[0][idx];// normal已初始化
    myfloat ny = edge.normal[1][idx];
    myfloat x_edge = edge.xy[0][idx];
    myfloat y_edge = edge.xy[1][idx];
    myfloat length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // 计算U_L，即edge U的左极限
    myfloat U_L[4]{};
    Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementL, U_L, para);
    // 计算ruvp_L
    myfloat ruvp_L[4]{};
    GPU::Math::U2ruvp(U_L, ruvp_L, para.constant.gamma);
    // 检查是否发散
    bool is_divergent = false;
    if (ruvp_L[3] < 0) {
        is_divergent = true;
        //std::cout << "p < 0, ";
    }
    if (ruvp_L[0] < 0) {
        is_divergent = true;
        //std::cout << "rho < 0, ";
    }
    if (ruvp_L[0] > 900) {
        is_divergent = true;
        //std::cout << "rho > 900, ";
    }
    if (is_divergent) {

    }


    //根据远场边界条件修正ruvp_L
    Space::Flux::modify_ruvpL_farField_kernel(nx, ny, ruvp_L, para);
    GPU::Math::ruvp2U(ruvp_L, U_L, para.constant.gamma);
    // 计算flux
    RiemannSolve_kernel(U_L, U_L, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::modify_ruvpL_farField_kernel(const myfloat nx, const myfloat ny, myfloat* ruvp, SDevicePara para) {
    // 该函数根据远场边界条件修正ruvp_L

    //边界元的Element_R为null，因此nxny一定朝外
    const myfloat rho = ruvp[0];
    const myfloat u = ruvp[1];
    const myfloat v = ruvp[2];
    const myfloat p = ruvp[3];
    const myfloat rho_inf = para.boundaryCondition_2D.ruvp_inf[0];
    const myfloat u_inf = para.boundaryCondition_2D.ruvp_inf[1];
    const myfloat v_inf = para.boundaryCondition_2D.ruvp_inf[2];
    const myfloat p_inf = para.boundaryCondition_2D.ruvp_inf[3];
    const myfloat gamma = para.constant.gamma;

    //坐标变换
    const myfloat rho_n = rho;
    const myfloat u_n = u * nx + v * ny;//u cost + v sint
    const myfloat v_n = -u * ny + v * nx;//此处有修改。法线向外为正 - u sint + v cost
    const myfloat p_n = p;
    const myfloat rho_n_inf = rho_inf;
    const myfloat p_n_inf = p_inf;
    const myfloat u_n_inf = u_inf * nx + v_inf * ny;
    const myfloat v_n_inf = -u_inf * ny + v_inf * nx;//此处有修改


    myfloat a_n2 = gamma * p_n / rho_n;
    if (a_n2 < 0) {
        printf("Error: a_n2 < 0 @GPU::Space::Flux::modify_ruvpL_farField_kernel\n");
        return;

    }
    const myfloat a_n = sqrt(a_n2);
    myfloat a_n_inf2 = gamma * p_n_inf / rho_n_inf;
    if (a_n_inf2 < 0) {
        printf("Error: a_n_inf2 < 0 @GPU::Space::Flux::modify_ruvpL_farField_kernel\n");
        return;

    }
    const myfloat a_n_inf = sqrt(a_n_inf2);

    myfloat v_n_tmp;
    myfloat u_n_tmp;
    myfloat rho_n_tmp;
    myfloat p_n_tmp;

    if (u_n_inf * u_n_inf < a_n_inf * a_n_inf) {//|Vn_inf|<|a_inf|，亚音速
        myfloat tmp_a = (gamma - 1.0) / 4.0 * (u_n - u_n_inf) + 0.5 * (a_n + a_n_inf);// 1/2 *( (ga-1)/2 * (u-uinf) + (a+ainf) )
        myfloat tmp_a2 = tmp_a * tmp_a;
        if (u_n_inf <= 0.0) {//亚音速入口，提3个边界条件
            myfloat tmp_s_inf = p_n_inf / pow(rho_n_inf, gamma);//根据边界rho和p求出熵
            rho_n_tmp = pow(tmp_a2 / tmp_s_inf / gamma, 1.0 / (gamma - 1.0));//边界
            u_n_tmp = 0.5 * (u_n + u_n_inf + 2.0 * a_n / (gamma - 1.0) - 2.0 * a_n_inf / (gamma - 1.0));
            v_n_tmp = v_n_inf;//边界
            p_n_tmp = tmp_a2 * rho_n_tmp / gamma;//边界
        }
        else {//亚音速出口，提1个边界条件 u_n_tmp
            myfloat tmp_s = p_n / pow(rho_n, gamma);//内 //pow(x,y):若x<0且y非整数，或者x=0且y<=0，将出现结果错误
            rho_n_tmp = pow(tmp_a2 / tmp_s / gamma, 1.0 / (gamma - 1.0));
            u_n_tmp = 0.5 * (u_n + u_n_inf + 2.0 * a_n / (gamma - 1.0) - 2.0 * a_n_inf / (gamma - 1.0));//边 //法向按公式计算u_n_i_m_η后，与u_n平均
            v_n_tmp = v_n;//内点
            p_n_tmp = tmp_a2 * rho_n_tmp / gamma;
        }
    }
    else {
        if (u_n_inf <= 0.0) {//超音速入口，提4个边界条件
            rho_n_tmp = rho_n_inf;
            u_n_tmp = u_n_inf;
            v_n_tmp = v_n_inf;
            p_n_tmp = p_n_inf;
        }
        else {//超音速出口，无需边界条件，所有关系由内点补充
            rho_n_tmp = rho_n;
            u_n_tmp = u_n;
            v_n_tmp = v_n;
            p_n_tmp = p_n;
        }
    }

    ruvp[0] = rho_n_tmp;
    ruvp[1] = u_n_tmp * nx - v_n_tmp * ny;//逆变换。这里用了自己的变换方式
    ruvp[2] = u_n_tmp * ny + v_n_tmp * nx;
    ruvp[3] = p_n_tmp;
}
__device__ void GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel(ElementSoA& element, ElementFieldSoA& elementField, EdgeSoA& edge, long idx, myfloat* flux, SDevicePara para) {
    //计算边中点坐标、坐标变换矩阵、坐标逆变换矩阵、边法线方向
    myfloat nx = edge.normal[0][idx];// normal已初始化
    myfloat ny = edge.normal[1][idx];
    myfloat x_edge = edge.xy[0][idx];
    myfloat y_edge = edge.xy[1][idx];
    myfloat length = edge.length[idx];
    long iElementL = edge.elementL[idx];
    long iElementR = edge.elementR[idx];


    // 计算edge坐标处的U值(分别用两侧单元的分布函数计算)，分为周期边界和内部边界两种情况
    myfloat U_L[4]{};
    Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementL, U_L, para);
    myfloat U_R[4]{};
    if (edge.setID[idx] != -1) {
        // 周期边界 应使用对称的edge的坐标计算U值
        //throw "Warning: periodic not implemented yet!\n";
        int ID = edge.ID[idx];
        int ID_pair = edge.periodicPair[ID];
        if (ID_pair < 0 || ID_pair >= edge.num_edge) {
            printf("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel\n");
            return;
            // 代办：异常处理。准备在GPUGlobalFunction中存一个全局变量ErrorCode，但是报错重定义
            // GPU::globalErrorCode = GPU::ErrorCode::periodicPairNotFound;
        }


        myfloat x_edge_pair = edge.xy[0][ID_pair];
        myfloat y_edge_pair = edge.xy[1][ID_pair];
        Flux::getUByXYandElementID_kernel(element, elementField, x_edge_pair, y_edge_pair, iElementR, U_R, para);

    }
    else {
        // 内部边界
        Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementR, U_R, para);
    }

    // 计算flux
    RiemannSolve_kernel(U_L, U_R, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::getUByXYandElementID_kernel(ElementSoA& element, ElementFieldSoA& elementField, myfloat x, myfloat y, int elementID, myfloat* U_dist, SDevicePara para) {
    // 该函数依据单元的重构函数，计算某坐标处的U值
    myfloat x_elementL = element.xy[0][elementID];
    myfloat y_elementL = element.xy[1][elementID];
    // 常量重构
    if (para.space.flag_reconstruct == _REC_constant) {
        for (int jValue = 0; jValue < 4; jValue++) {
            U_dist[jValue] = elementField.U[jValue][elementID];
        }
    }
    // 线性重构 根据梯度Ux、Uy计算点(xpoint,ypoint)处_U值
    else if (para.space.flag_reconstruct == _REC_linear) {
        for (int jValue = 0; jValue < 4; jValue++) {
            myfloat& U_elementL = elementField.U[jValue][elementID];
            myfloat& Ux_elementL = elementField.Ux[jValue][elementID];
            myfloat& Uy_elementL = elementField.Uy[jValue][elementID];
            U_dist[jValue] = U_elementL + Ux_elementL * (x - x_elementL) + Uy_elementL * (y - y_elementL);
        }
        // 若数据异常，则常量重构
        myfloat ruvp[4]{};
        GPU::Math::U2ruvp(U_dist, ruvp, para.constant.gamma);
        if (Space::Restrict::outOfRange(ruvp)) {
            for (int jValue = 0; jValue < 4; jValue++) {
                U_dist[jValue] = elementField.U[jValue][elementID];
            }
        }
    }
    else {
        printf("Error: invalid reconstruct type.\n");
    }
}

__device__ void GPU::Space::Flux::RiemannSolve_kernel(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny,
    const myfloat length, myfloat* flux, SDevicePara para) {
    myfloat faceNormal[2]{ nx,ny };

    switch (para.inviscidFluxMethod.flux_conservation_scheme) {
    case _SOL_LocalLaxFriedrichs:
        GPU::Space::Convection::LocalLaxFriedrichs2d(UL, UR, nx, ny, length, flux, para.constant.gamma);
        break;
    case _SOL_Roe:
        GPU::Space::Convection::ConvectRoeCommon2d(UL, UR, faceNormal, length, flux, para);
        break;
    default:
        break;
    }
}


__global__ void GPU::Space::Flux::calculateFluxKernel(ElementSoA element, ElementFieldSoA elementField, EdgeSoA edge, BoundarySetMap boundary_device, SDevicePara para) {
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
        GPU::Space::Flux::getEdgeFlux_wallNonViscous_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //无粘固壁
    case _BC_wall_nonViscous:
        GPU::Space::Flux::getEdgeFlux_wallNonViscous_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //无滑移绝热
    case _BC_wall_adiabat:
        // TODO: 实现无滑移绝热边界条件的计算
        GPU::Space::Flux::getEdgeFlux_wall_adiabat_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //入口
    case _BC_inlet:
        GPU::Space::Flux::getEdgeFlux_farField_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //出口
    case _BC_outlet:
        GPU::Space::Flux::getEdgeFlux_farField_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //远场
    case _BC_inf:
        GPU::Space::Flux::getEdgeFlux_farField_kernel(element, elementField, edge, iedge, flux, para);
        break;
    default:// 内部：bType=-1
        if (elementR != -1) {
            // 周期和内部 统一处理
            GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel(element, elementField, edge, iedge, flux, para);
        }
    }

    // 更新单元数值通量 由于是多线程运行，求和时会发生数据竞争
    for (int j = 0; j < 4; j++) {
        elementField.Flux[j][elementL] += flux[j];
        // 内部边界
        if (bType == -1) {
            // 此处不是elementR==-1，防止周期边界算2次
            elementField.Flux[j][elementR] -= flux[j];
        }
    }
}


void GPU::Space::Flux::calculateFlux(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, SDevicePara& infPara_device) {
    //LogWriter::logAndPrintError("Not completed. @GPU::Space::Flux::calculateFlux\n");
    //exit(-1);// 未完成 getEdgeFlux_inner_and_periodic
    int block_size = 512;
    int grid_size = (edge_device.num_edge + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    calculateFluxKernel <<<grid, block >>> (element_device, elementField_device, edge_device, boundary_device, infPara_device);
}


void GPU::Space::Flux::calculateFluxDevice(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, SDevicePara& infPara_device) {
    // TODO: Implement this function.
    resetElementFlux_1(elementField_device);
    cudaDeviceSynchronize();
    catchCudaErrorAndExit();
    calculateFlux(element_device, elementField_device, edge_device, boundary_device, infPara_device);
}


// 以下为新代码


__global__ void resetElementFluxKernel_2(GPU::ElementFieldSoA elementField_device) {
    // --- 核函数 --- 

    // 获取id，判断是否有效
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& i = id;
    if (i >= elementField_device.num || i < 0) return;

    // 给单元数值通量赋值
    for (int jValue = 0; jValue < 4; jValue++) {
        elementField_device.Flux[jValue][i] = 0;//单元数值通量Flux清零，为后面加减做准备 Invalid __global__ write of size 8 bytes
        // f->elements[ie].deltaeig = 0;//每一轮deltaeig清零 针对Roe格式 目前不需要
        // 梯度已经计算，因此不需要再次计算
    }
}

__global__ void setVectorToZero(integer length, myfloat* vector) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= length)return;

}

void resetElementFlux_2(GPU::ElementFieldSoA& elementField_device) {
    // 初始化单元数值通量
    int block_size = 512;// 最好是128 256 512
    int grid_size = (elementField_device.num + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    resetElementFluxKernel_2 << <grid, block >> > (elementField_device);
}


void GPU::Space::Flux::calculateFluxDevice_2(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, SDevicePara& para) {
    LogWriter::logAndPrintError("unimplemented @ GPU::Space::Flux::calculateFluxDevice_2");
    exit(-1);

    /*
    https://forums.developer.nvidia.com/t/synchronization-between-kernel-calls/23336
    两次核函数调用之间无需加cudaDeviceSynchronize()，因为 Kernels in the same stream are invoked sequentially, so you don’t need any extra synchronization in this example.
    If your kernel returns unexpected results, it must be for some other reason.
    */

    integer num = elementField_device.num;
    // 通量清零
    for (int i = 0; i < 4; i++) {
        cudaMemset(elementField_device.Flux[i], 0, num * sizeof(myfloat));
    }
    catchCudaErrorAndExit();


}
