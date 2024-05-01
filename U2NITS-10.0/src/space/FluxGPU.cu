#include "FluxGPU.h"
#include "../FVM_2D.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../math/MathGPU.h"
#include "convection/ConvectionGPU.h"
#include "../output/LogWriter.h"
#include "restrict/RestrictGPU.h"

__global__ void resetElementFluxKernel_1(GPU::ElementFieldSoA elementField_device) {
    // --- �˺��� --- 

    // ��ȡid���ж��Ƿ���Ч
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& i = id;
    if (i >= elementField_device.num || i < 0) return;

    // ����Ԫ��ֵͨ����ֵ
    for (int jValue = 0; jValue < 4; jValue++) {
        elementField_device.Flux[jValue][i] = 0;//��Ԫ��ֵͨ��Flux���㣬Ϊ����Ӽ���׼�� Invalid __global__ write of size 8 bytes
        // f->elements[ie].deltaeig = 0;//ÿһ��deltaeig���� ���Roe��ʽ Ŀǰ����Ҫ
        // �ݶ��Ѿ����㣬��˲���Ҫ�ٴμ���
    }
}

void resetElementFlux_1(GPU::ElementFieldSoA& elementField_device) {
    // ��ʼ����Ԫ��ֵͨ��
    int block_size = 512;// �����128 256 512
    int grid_size = (elementField_device.num + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    resetElementFluxKernel_1 <<<grid, block>>> (elementField_device);
}

__device__ void GPU::Space::Flux::getEdgeFlux_wallNonViscous_kernel(ElementSoA& element, ElementFieldSoA& elementField, EdgeSoA& edge, long idx, myfloat* flux, SDevicePara para) {
    //δ�����;
    // ������idx�൱��iedge
    // �ú����ֱ����ߵ���������U�ļ��ޣ�Ȼ�����flux

    myfloat nx = edge.normal[0][idx];// normal�ѳ�ʼ��
    myfloat ny = edge.normal[1][idx];
    myfloat x_edge = edge.xy[0][idx];
    myfloat y_edge = edge.xy[1][idx];
    myfloat length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // ����U_L����edge U������
    myfloat U_L[4]{};
    Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementL, U_L, para);

    // ����U_R �öԳ���
    myfloat U_R[4]{};
    myfloat uL = U_L[1] / U_L[0];
    myfloat vL = U_L[2] / U_L[0];
    myfloat uR = uL - 2 * (uL * nx + vL * ny) * nx;
    myfloat vR = vL - 2 * (uL * nx + vL * ny) * ny;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    // ����flux
    RiemannSolve_kernel(U_L, U_R, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::getEdgeFlux_wall_adiabat_kernel(ElementSoA& element, ElementFieldSoA& elementField, EdgeSoA& edge, long idx, myfloat* flux, SDevicePara para) {
    // �̱ڱ߽硣������idx�൱��iedge

    myfloat nx = edge.normal[0][idx];// normal�ѳ�ʼ��
    myfloat ny = edge.normal[1][idx];
    myfloat x_edge = edge.xy[0][idx];
    myfloat y_edge = edge.xy[1][idx];
    myfloat length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // ����U_L����edge U������
    myfloat U_L[4]{};
    Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementL, U_L, para);

    // ����U_R �öԳ���
    myfloat U_R[4]{};
    myfloat uL = U_L[1] / U_L[0];
    myfloat vL = U_L[2] / U_L[0];
    // ��������ϵ�� UnL = uL * nx + vL * ny��VnL = -uL * ny + vL * nx
    // UnR = - UnL��VnR = -VnL
    // ԭ����ϵ�� uR = unR*nx-vnR*ny�� vR=unR*ny+vnR*nx
    // ��� uR = -uL ������ֱ��ȡ��

    // ��UNITS��BC_WALL.f90���㷨���ٶ�ȡ�������ǽ����Ȼ����
    // ���� �����ͼ��û������ճ��
    //myfloat uR = uL - 2 * (uL * nx + vL * ny) * nx;
    //myfloat vR = vL - 2 * (uL * nx + vL * ny) * ny;
    myfloat uR = -uL;
    myfloat vR = -vL;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    // ����flux
    RiemannSolve_kernel(U_L, U_R, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::getEdgeFlux_farField_kernel(ElementSoA& element, ElementFieldSoA& elementField, EdgeSoA& edge, long idx, myfloat* flux, SDevicePara para) {
    // �������ܣ�����Զ���߽���������߽���ֵͨ�� �����

    myfloat nx = edge.normal[0][idx];// normal�ѳ�ʼ��
    myfloat ny = edge.normal[1][idx];
    myfloat x_edge = edge.xy[0][idx];
    myfloat y_edge = edge.xy[1][idx];
    myfloat length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // ����U_L����edge U������
    myfloat U_L[4]{};
    Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementL, U_L, para);
    // ����ruvp_L
    myfloat ruvp_L[4]{};
    GPU::Math::U2ruvp(U_L, ruvp_L, para.constant.gamma);
    // ����Ƿ�ɢ
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


    //����Զ���߽���������ruvp_L
    Space::Flux::modify_ruvpL_farField_kernel(nx, ny, ruvp_L, para);
    GPU::Math::ruvp2U(ruvp_L, U_L, para.constant.gamma);
    // ����flux
    RiemannSolve_kernel(U_L, U_L, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::modify_ruvpL_farField_kernel(const myfloat nx, const myfloat ny, myfloat* ruvp, SDevicePara para) {
    // �ú�������Զ���߽���������ruvp_L

    //�߽�Ԫ��Element_RΪnull�����nxnyһ������
    const myfloat rho = ruvp[0];
    const myfloat u = ruvp[1];
    const myfloat v = ruvp[2];
    const myfloat p = ruvp[3];
    const myfloat rho_inf = para.boundaryCondition_2D.ruvp_inf[0];
    const myfloat u_inf = para.boundaryCondition_2D.ruvp_inf[1];
    const myfloat v_inf = para.boundaryCondition_2D.ruvp_inf[2];
    const myfloat p_inf = para.boundaryCondition_2D.ruvp_inf[3];
    const myfloat gamma = para.constant.gamma;

    //����任
    const myfloat rho_n = rho;
    const myfloat u_n = u * nx + v * ny;//u cost + v sint
    const myfloat v_n = -u * ny + v * nx;//�˴����޸ġ���������Ϊ�� - u sint + v cost
    const myfloat p_n = p;
    const myfloat rho_n_inf = rho_inf;
    const myfloat p_n_inf = p_inf;
    const myfloat u_n_inf = u_inf * nx + v_inf * ny;
    const myfloat v_n_inf = -u_inf * ny + v_inf * nx;//�˴����޸�


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

    if (u_n_inf * u_n_inf < a_n_inf * a_n_inf) {//|Vn_inf|<|a_inf|��������
        myfloat tmp_a = (gamma - 1.0) / 4.0 * (u_n - u_n_inf) + 0.5 * (a_n + a_n_inf);// 1/2 *( (ga-1)/2 * (u-uinf) + (a+ainf) )
        myfloat tmp_a2 = tmp_a * tmp_a;
        if (u_n_inf <= 0.0) {//��������ڣ���3���߽�����
            myfloat tmp_s_inf = p_n_inf / pow(rho_n_inf, gamma);//���ݱ߽�rho��p�����
            rho_n_tmp = pow(tmp_a2 / tmp_s_inf / gamma, 1.0 / (gamma - 1.0));//�߽�
            u_n_tmp = 0.5 * (u_n + u_n_inf + 2.0 * a_n / (gamma - 1.0) - 2.0 * a_n_inf / (gamma - 1.0));
            v_n_tmp = v_n_inf;//�߽�
            p_n_tmp = tmp_a2 * rho_n_tmp / gamma;//�߽�
        }
        else {//�����ٳ��ڣ���1���߽����� u_n_tmp
            myfloat tmp_s = p_n / pow(rho_n, gamma);//�� //pow(x,y):��x<0��y������������x=0��y<=0�������ֽ������
            rho_n_tmp = pow(tmp_a2 / tmp_s / gamma, 1.0 / (gamma - 1.0));
            u_n_tmp = 0.5 * (u_n + u_n_inf + 2.0 * a_n / (gamma - 1.0) - 2.0 * a_n_inf / (gamma - 1.0));//�� //���򰴹�ʽ����u_n_i_m_�Ǻ���u_nƽ��
            v_n_tmp = v_n;//�ڵ�
            p_n_tmp = tmp_a2 * rho_n_tmp / gamma;
        }
    }
    else {
        if (u_n_inf <= 0.0) {//��������ڣ���4���߽�����
            rho_n_tmp = rho_n_inf;
            u_n_tmp = u_n_inf;
            v_n_tmp = v_n_inf;
            p_n_tmp = p_n_inf;
        }
        else {//�����ٳ��ڣ�����߽����������й�ϵ���ڵ㲹��
            rho_n_tmp = rho_n;
            u_n_tmp = u_n;
            v_n_tmp = v_n;
            p_n_tmp = p_n;
        }
    }

    ruvp[0] = rho_n_tmp;
    ruvp[1] = u_n_tmp * nx - v_n_tmp * ny;//��任�����������Լ��ı任��ʽ
    ruvp[2] = u_n_tmp * ny + v_n_tmp * nx;
    ruvp[3] = p_n_tmp;
}
__device__ void GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel(ElementSoA& element, ElementFieldSoA& elementField, EdgeSoA& edge, long idx, myfloat* flux, SDevicePara para) {
    //������е����ꡢ����任����������任���󡢱߷��߷���
    myfloat nx = edge.normal[0][idx];// normal�ѳ�ʼ��
    myfloat ny = edge.normal[1][idx];
    myfloat x_edge = edge.xy[0][idx];
    myfloat y_edge = edge.xy[1][idx];
    myfloat length = edge.length[idx];
    long iElementL = edge.elementL[idx];
    long iElementR = edge.elementR[idx];


    // ����edge���괦��Uֵ(�ֱ������൥Ԫ�ķֲ���������)����Ϊ���ڱ߽���ڲ��߽��������
    myfloat U_L[4]{};
    Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementL, U_L, para);
    myfloat U_R[4]{};
    if (edge.setID[idx] != -1) {
        // ���ڱ߽� Ӧʹ�öԳƵ�edge���������Uֵ
        //throw "Warning: periodic not implemented yet!\n";
        int ID = edge.ID[idx];
        int ID_pair = edge.periodicPair[ID];
        if (ID_pair < 0 || ID_pair >= edge.num_edge) {
            printf("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel\n");
            return;
            // ���죺�쳣����׼����GPUGlobalFunction�д�һ��ȫ�ֱ���ErrorCode�����Ǳ����ض���
            // GPU::globalErrorCode = GPU::ErrorCode::periodicPairNotFound;
        }


        myfloat x_edge_pair = edge.xy[0][ID_pair];
        myfloat y_edge_pair = edge.xy[1][ID_pair];
        Flux::getUByXYandElementID_kernel(element, elementField, x_edge_pair, y_edge_pair, iElementR, U_R, para);

    }
    else {
        // �ڲ��߽�
        Flux::getUByXYandElementID_kernel(element, elementField, x_edge, y_edge, iElementR, U_R, para);
    }

    // ����flux
    RiemannSolve_kernel(U_L, U_R, nx, ny, length, flux, para);
}
__device__ void GPU::Space::Flux::getUByXYandElementID_kernel(ElementSoA& element, ElementFieldSoA& elementField, myfloat x, myfloat y, int elementID, myfloat* U_dist, SDevicePara para) {
    // �ú������ݵ�Ԫ���ع�����������ĳ���괦��Uֵ
    myfloat x_elementL = element.xy[0][elementID];
    myfloat y_elementL = element.xy[1][elementID];
    // �����ع�
    if (para.space.flag_reconstruct == _REC_constant) {
        for (int jValue = 0; jValue < 4; jValue++) {
            U_dist[jValue] = elementField.U[jValue][elementID];
        }
    }
    // �����ع� �����ݶ�Ux��Uy�����(xpoint,ypoint)��_Uֵ
    else if (para.space.flag_reconstruct == _REC_linear) {
        for (int jValue = 0; jValue < 4; jValue++) {
            myfloat& U_elementL = elementField.U[jValue][elementID];
            myfloat& Ux_elementL = elementField.Ux[jValue][elementID];
            myfloat& Uy_elementL = elementField.Uy[jValue][elementID];
            U_dist[jValue] = U_elementL + Ux_elementL * (x - x_elementL) + Uy_elementL * (y - y_elementL);
        }
        // �������쳣�������ع�
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
    // δ���

    // ���ڱ߽磬EdgeSoA��periodicPair��Ա��������ɣ��Ѿ�������edge_device�У������ٴ���
    // ��ȡ�߽�����ʱ������ʹ��boundary_host

    // ��ȡid���ж��Ƿ���Ч
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& iedge = id;
    if (iedge >= edge.num_edge || iedge < 0) return;

    // ÿ���߼�����ճ��ֵͨ��

    //FVM_2D* f = FVM_2D::getInstance();
    int elementL = edge.elementL[iedge];
    int elementR = edge.elementR[iedge];
    int setID = edge.setID[iedge];// �ڲ�edge�������κ�set�����setIDΪ-1 setID�ĳ�ʼ����readContinueFile
    //int boundaryNum = f->boundaryManager.boundaries.size();
    int boundaryNum = boundary_device.size;
    int bType = -1;
    if (setID != -1) {
        bType = boundary_device.type[setID - 1];
    }
    double flux[4]{};

    switch (bType) {
        //�Գơ���ŷ�����̣��൱����ճ�̱�
    case _BC_symmetry:
        GPU::Space::Flux::getEdgeFlux_wallNonViscous_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //��ճ�̱�
    case _BC_wall_nonViscous:
        GPU::Space::Flux::getEdgeFlux_wallNonViscous_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //�޻��ƾ���
    case _BC_wall_adiabat:
        // TODO: ʵ���޻��ƾ��ȱ߽������ļ���
        GPU::Space::Flux::getEdgeFlux_wall_adiabat_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //���
    case _BC_inlet:
        GPU::Space::Flux::getEdgeFlux_farField_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //����
    case _BC_outlet:
        GPU::Space::Flux::getEdgeFlux_farField_kernel(element, elementField, edge, iedge, flux, para);
        break;
        //Զ��
    case _BC_inf:
        GPU::Space::Flux::getEdgeFlux_farField_kernel(element, elementField, edge, iedge, flux, para);
        break;
    default:// �ڲ���bType=-1
        if (elementR != -1) {
            // ���ں��ڲ� ͳһ����
            GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel(element, elementField, edge, iedge, flux, para);
        }
    }

    // ���µ�Ԫ��ֵͨ�� �����Ƕ��߳����У����ʱ�ᷢ�����ݾ���
    for (int j = 0; j < 4; j++) {
        elementField.Flux[j][elementL] += flux[j];
        // �ڲ��߽�
        if (bType == -1) {
            // �˴�����elementR==-1����ֹ���ڱ߽���2��
            elementField.Flux[j][elementR] -= flux[j];
        }
    }
}


void GPU::Space::Flux::calculateFlux(ElementSoA& element_device, ElementFieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, SDevicePara& infPara_device) {
    //LogWriter::logAndPrintError("Not completed. @GPU::Space::Flux::calculateFlux\n");
    //exit(-1);// δ��� getEdgeFlux_inner_and_periodic
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


// ����Ϊ�´���


__global__ void resetElementFluxKernel_2(GPU::ElementFieldSoA elementField_device) {
    // --- �˺��� --- 

    // ��ȡid���ж��Ƿ���Ч
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& i = id;
    if (i >= elementField_device.num || i < 0) return;

    // ����Ԫ��ֵͨ����ֵ
    for (int jValue = 0; jValue < 4; jValue++) {
        elementField_device.Flux[jValue][i] = 0;//��Ԫ��ֵͨ��Flux���㣬Ϊ����Ӽ���׼�� Invalid __global__ write of size 8 bytes
        // f->elements[ie].deltaeig = 0;//ÿһ��deltaeig���� ���Roe��ʽ Ŀǰ����Ҫ
        // �ݶ��Ѿ����㣬��˲���Ҫ�ٴμ���
    }
}

__global__ void setVectorToZero(integer length, myfloat* vector) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= length)return;

}

void resetElementFlux_2(GPU::ElementFieldSoA& elementField_device) {
    // ��ʼ����Ԫ��ֵͨ��
    int block_size = 512;// �����128 256 512
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
    ���κ˺�������֮�������cudaDeviceSynchronize()����Ϊ Kernels in the same stream are invoked sequentially, so you don��t need any extra synchronization in this example.
    If your kernel returns unexpected results, it must be for some other reason.
    */

    integer num = elementField_device.num;
    // ͨ������
    for (int i = 0; i < 4; i++) {
        cudaMemset(elementField_device.Flux[i], 0, num * sizeof(myfloat));
    }
    catchCudaErrorAndExit();


}
