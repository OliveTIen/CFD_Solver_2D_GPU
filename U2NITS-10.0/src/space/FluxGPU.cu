#include "FluxGPU.h"
#include "../FVM_2D.h"
#include "../gpu/GPUGlobalFunction.h"
#include "../math/MathGPU.h"
#include "convection/ConvectionGPU.h"
#include "../output/LogWriter.h"
#include "restrict/RestrictGPU.h"
#include "../solvers/SolverDataGetter.h"
#include "../solvers/GPUSolver2.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"
#include "viscous_flux/ViscousFluxGPU.h"


// ��Ԫͨ������
inline void clear_element_flux_device() {
    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementFieldSoA& elementField_device = solver->elementField_device;

    integer num = elementField_device.num;
    for (int i = 0; i < 4; i++) {
        cudaMemset(elementField_device.Flux[i], 0, num * sizeof(myfloat));
    }
    
    getLastCudaError("clear_element_flux_device failed.");
}

__device__ inline bool isUpStreamOfShock_atBoundary_static_inline(double x, double y, double _shock_x, double _shock_y, double _shock_speed, double _t_RK, const double _sqrt3) {
    // CBoundaryDoubleShockReflect�ĳ�Ա�����ľ�̬�汾
    double right = _shock_x + (y + _shock_speed * _t_RK) / _sqrt3;
    if (x < right) return true;
    else return false;
}

__device__ inline void get_U_reconstruct_singleValue(myfloat x, myfloat y, myfloat& U_dist, myint i_e, const myfloat* x_e, const myfloat* y_e, const myfloat* U_e, const myfloat* Ux_e, const myfloat* Uy_e, int flag_reconstruct) {
    if (flag_reconstruct == _REC_constant) {
        U_dist = U_e[i_e];
    }
    else if (flag_reconstruct == _REC_linear) {
        U_dist = U_e[i_e] + Ux_e[i_e] * (x - x_e[i_e]) + Uy_e[i_e] * (y - y_e[i_e]);
    }
    else {
        printf("invalid flag_reconstruct. @get_U_reconstruct_singleValue.\n");
        U_dist = U_e[i_e];
    }
}

__device__ void get_Uvector_reconstruct_2_device(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, myfloat x, myfloat y, myint iElement, myfloat* U_dist, int flag_reconstruct, myfloat gamma) {
    // ���ݵ�Ԫ�ֲ��������õ�ĳ��Uֵ
    for (int i = 0; i < 4; i++) {
        get_U_reconstruct_singleValue(x, y, U_dist[i], iElement, element.xy[0], element.xy[1], elementField.U[i], elementField.Ux[i], elementField.Uy[i], flag_reconstruct);
    }

    // �������쳣�������ع�
    myfloat ruvp[4]{};
    GPU::Math::U2ruvp(U_dist, ruvp, gamma);
    if (GPU::Space::Restrict::outOfRange(ruvp)) {
        for (int i = 0; i < 4; i++) {
            U_dist[i] = elementField.U[i][iElement];
        }
    }
}

__device__ void get_UR_wallNonViscous_device(const myfloat* U_L, myfloat* U_R, myfloat nx, myfloat ny) {
    // ���Ʊ��棬�����ٶ�ȡ���������ٶ����
    myfloat uxL = U_L[1] / U_L[0];
    myfloat uyL = U_L[2] / U_L[0];
    myfloat unL = uxL * nx + uyL * ny;// xyΪȫ������ϵ��tnΪ����-��������ϵ
    myfloat utL = uxL * ny - uyL * nx;
    myfloat unR = -unL;
    myfloat utR = utL;
    myfloat uxR = unR * nx + utR * ny;
    myfloat uyR = unR * ny - utR * nx;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uxR;
    U_R[2] = U_R[0] * uyR;
    U_R[3] = U_L[3];
}

__device__ void get_UR_wall_adiabat_device(const myfloat* U_L, myfloat* U_R, myfloat nx, myfloat ny) {
    // �޻��ƾ��ȱ��棬���������ٶȶ�ȡ��
    U_R[0] = U_L[0];
    U_R[1] = -U_R[1];
    U_R[2] = -U_R[2];
    U_R[3] = U_L[3];
}

__device__ void get_UR_farField_device(const myfloat* U_L, const myfloat* ruvp_inf, myfloat* U_R, myfloat nx, myfloat ny, myfloat gamma) {
    /*
    ����UNITs���������μ�
    �����������Ϊ���������ڵ����Ʊ߽��Ĺ���ֵ��
    �ٽ����߽����ֵ�������ֵ���������ݹ�ϵ
    ����߽����٣������õ�ѹ����ֵ��
    */
    myfloat ga1 = gamma - 1;
    myfloat rho_f = ruvp_inf[0];
    myfloat u_f = ruvp_inf[1];
    myfloat v_f = ruvp_inf[2];
    myfloat p_f = ruvp_inf[3];

    myfloat ruvp_e[4]{};
    GPU::Math::U2ruvp(U_L, ruvp_e, gamma);
    myfloat rho_e = ruvp_e[0];
    myfloat u_e = ruvp_e[1];
    myfloat v_e = ruvp_e[2];
    myfloat p_e = ruvp_e[3];

    // �������ݹ�ϵ
    myfloat a2_f = gamma * p_f / rho_f;
    myfloat af = sqrt(a2_f);// Զ������
    myfloat Ma2_f = (u_f * u_f + v_f * v_f) / a2_f;// Զ�������ƽ��
    myfloat ae = sqrt(gamma * p_e / rho_e);// �ڵ�����
    myfloat qnf = u_f * nx + v_f * ny;// Զ�������ٶ�
    myfloat qne = u_e * nx + v_e * ny;// �ڵ㷨���ٶ�
    myfloat rf = qnf - 2.0 * af / ga1;// ��һ������ R2 = u0n - 2*a0n/(gamma-1)
    myfloat re = qne + 2.0 * ae / ga1;// �ڶ������� u0n - 2*a0n/(gamma-1)
    myfloat qn = 0.5 * (re + rf);
    myfloat as = ga1 * (re - rf) / 4.0;
    myfloat dnt = 0;// �����˶��ٶ�

    // �о� q<=-as����������ڣ�-as<q<0����������ڣ�0<=q<as�������ٳ��ڣ�as<=q�������ٳ���
    myfloat rho0 = 1, u0 = 0, v0 = 0, p0 = 1;
    if (qn <= -as) {
        // ���������
        rho0 = rho_f;
        u0 = u_f;
        v0 = v_f;
        p0 = p_f;
    }
    else if (qn < 0) {
        // -as < qn < 0 ע��C++�в�������дС�ں�
        // ��������� ����UNITs��
        myfloat son = p_f / pow(rho_f, gamma);// ��
        myfloat qtx = u_f - qnf * nx;
        myfloat qty = v_f - qnf * ny;

        rho0 = pow((as * as / son / gamma), 1.0 / ga1);// ����S=p/rho^gamma��rho=(p/S)^(1/gamma)=(a2/gamma/S)^(1/gamma)
        u0 = qtx + qn * nx;
        v0 = qty + qn * ny;
        p0 = as * as * rho0 / gamma;// p=rho*R*t=rho/gamma*gamma*R*t=rho/gamma*a2
    }
    else if (qn < as) {
        // 0 <= qn < as
        // �����ٳ��� 
        myfloat son = p_e / pow(rho_e, gamma);
        myfloat qtx = u_e - qne * nx;
        myfloat qty = v_e - qne * ny;

        // ����������������ͬ
        rho0 = pow((as * as / son / gamma), 1.0 / ga1);
        u0 = qtx + qn * nx;
        v0 = qty + qn * ny;
        p0 = as * as * rho0 / gamma;
    }
    else {
        // as <= q
        // �����ٳ���
        rho0 = rho_e;
        u0 = u_e;
        v0 = v_e;
        p0 = p_e;
    }
    myfloat ruvp0[4]{ rho0,u0,v0,p0 };
    // ����flux
    if (GPU::Space::Restrict::outOfRange(ruvp0)) {
        printf("Out of range @get_UR_farField \n");
    }
    GPU::Math::ruvp2U(ruvp0, U_R, gamma);
}

__device__ void get_UR_inner_and_periodic_device(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, GPU::EdgeSoA& edge, myint iEdge, myint iElementR, myfloat x_edge, myfloat y_edge, myfloat* U_R, int inviscid_flux_method_flag_reconstruct, myfloat gamma) {
    if (edge.setID[iEdge] != -1) {
        // ���ڱ߽� Ӧʹ�öԳƵ�edge���������Uֵ
        int ID = edge.ID[iEdge];
        int ID_pair = edge.periodicPair[ID];
        if (ID_pair < 0 || ID_pair >= edge.num_edge) {
            printf("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel\n");
            return;
        }

        myfloat x_edge_pair = edge.xy[0][ID_pair];
        myfloat y_edge_pair = edge.xy[1][ID_pair];
        get_Uvector_reconstruct_2_device(element, elementField, x_edge_pair, y_edge_pair, iElementR, U_R, inviscid_flux_method_flag_reconstruct, gamma);
    }
    else {
        // �ڲ��߽�
        get_Uvector_reconstruct_2_device(element, elementField, x_edge, y_edge, iElementR, U_R, inviscid_flux_method_flag_reconstruct, gamma);
    }

}

__device__ void get_flux_RiemannSolve_device(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny, const myfloat length, myfloat* flux, const int flux_conservation_scheme, myfloat gamma, myfloat rcpcv) {
    // ���������������UL��UR���flux
    myfloat faceNormal[2]{ nx,ny };

    switch (flux_conservation_scheme) {// GlobalPara::inviscid_flux_method::flux_conservation_scheme
    case _SOL_LocalLaxFriedrichs:
        GPU::Space::Convection::LocalLaxFriedrichs2d(UL, UR, nx, ny, length, flux, gamma);
        break;
    case _SOL_Roe:
        GPU::Space::Convection::ConvectRoeCommon2d(UL, UR, faceNormal, length, flux, gamma, rcpcv);
        break;
    default:
        printf("invalid scheme @get_flux_RiemannSolve");
        break;
    }
}

__global__ void edge_convection_flux_device_kernel(GPU::ElementSoA element_device, GPU::ElementFieldSoA elementField_device, GPU::EdgeSoA edge_device, GPU::EdgeFieldSoA edgeField_device, GPU::BoundarySetMap boundary_device,
    int flux_conservation_scheme, int inviscid_flux_method_flag_reconstruct, myfloat* inf_ruvp, myfloat* inlet_ruvp, myfloat* outlet_ruvp, myfloat gamma, myfloat rcpcv,
    myfloat dsr_shock_x, myfloat dsr_shock_y, myfloat dsr_shock_speed, myfloat dsr_t_RK, const myfloat _sqrt3) {
    

    const int iEdge = blockIdx.x * blockDim.x + threadIdx.x;
    if (iEdge >= edge_device.num_edge || iEdge < 0) return;

    int setID = edge_device.setID[iEdge];// �ڲ�edge�������κ�set�����setIDΪ-1 setID�ĳ�ʼ����readContinueFile
    int bType = -1;
    if (setID != -1) {
        // ��iedge��edge�ı߽�����
        bType = boundary_device.type[setID - 1];
    }

    myfloat nx = edge_device.normal[0][iEdge];
    myfloat ny = edge_device.normal[1][iEdge];
    myfloat x_edge = edge_device.xy[0][iEdge];
    myfloat y_edge = edge_device.xy[1][iEdge];
    myfloat length = edge_device.length[iEdge];
    myint iElementL = edge_device.elementL[iEdge];
    myint iElementR = edge_device.elementR[iEdge];
    bool isUpStream_doubleShockReflect = isUpStreamOfShock_atBoundary_static_inline(x_edge, y_edge, dsr_shock_x, dsr_shock_y, dsr_shock_speed, dsr_t_RK, _sqrt3);

    myfloat U_L[4]{};
    myfloat U_R[4]{};
    myfloat flux[4]{};
    get_Uvector_reconstruct_2_device(element_device, elementField_device, x_edge, y_edge, iElementL, U_L, inviscid_flux_method_flag_reconstruct, gamma);// get U_L

    switch (bType) {
    case _BC_symmetry:
        // �Գơ���ŷ�����̣��൱����ճ�̱�
        get_UR_wallNonViscous_device(U_L, U_R, nx, ny);
        break;
    case _BC_wall_nonViscous:
        // ��ճ�̱�
        get_UR_wallNonViscous_device(U_L, U_R, nx, ny);
        break;
    case _BC_wall_adiabat:
        // �޻��ƾ���
        get_UR_wall_adiabat_device(U_L, U_R, nx, ny);
        break;
    case _BC_inlet:
        // ���
        get_UR_farField_device(U_L, inlet_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_outlet:
        // ����
        get_UR_farField_device(U_L, outlet_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_inf:
        // Զ��
        get_UR_farField_device(U_L, inf_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_doubleShockReflect:
        // ˫��շ���
        if (isUpStream_doubleShockReflect) {
            get_UR_farField_device(U_L, inlet_ruvp, U_R, nx, ny, gamma);
        }
        else {
            get_UR_farField_device(U_L, outlet_ruvp, U_R, nx, ny, gamma);
        }
        break;

    default:// �ڲ���bType=-1���߽磺bTypeȡ_BC_periodic_0��_BC_periodic_9����6100-6109
        if (iElementR != -1) {
            // ���ں��ڲ� ͳһ����
            get_UR_inner_and_periodic_device(element_device, elementField_device, edge_device, iEdge, iElementR, x_edge, y_edge, U_R, inviscid_flux_method_flag_reconstruct, gamma);
        }
    }

    get_flux_RiemannSolve_device(U_L, U_R, nx, ny, length, flux, flux_conservation_scheme, gamma, rcpcv);

    for (int j = 0; j < 4; j++) {
        edgeField_device.Flux[j][iEdge] = flux[j];
    }
}

__global__ void add_edge_flux_to_element_flux_kernel(GPU::ElementSoA element_device, GPU::ElementFieldSoA elementField_device, GPU::EdgeSoA edge_device, GPU::EdgeFieldSoA edgeField_device) {
    
    const int iElement = blockIdx.x * blockDim.x + threadIdx.x;
    if (iElement >= element_device.num_element || iElement < 0) return;

    myfloat volumeC = element_device.volume[iElement];// ���
    // ��ʼ��edgeID �������鷶Χ�ķǷ�ֵ��Ϊ-1
    int edgeID[4]{ -1,-1,-1,-1 };// ��iElement��element�ĵ�i���ߵ�ID
    for (int i = 0; i < 4; i++) {
        edgeID[i] = element_device.edges[i][iElement];
        if (edgeID[i] < 0 || edgeID[i] >= edge_device.num_edge) { // has: iEdge >= 0 && iEdge < num_edge
            edgeID[i] = -1;
        }
    }
    // ��ʼ���߳���edgeSign
    int edgeSign[4]{ 0,0,0,0 };// 1��ʾ�߳��⣬-1��ʾ�߳���
    for (int i = 0; i < 4; i++) {
        int edgeIDi = edgeID[i];
        if (edgeIDi == -1) {
            continue;
        }
        if (iElement != edge_device.elementR[edgeIDi]) {
            edgeSign[i] = 1;// currentElement=elementL������
        }
        else {
            edgeSign[i] = -1;// currentElement=elementR������
        }
    }
    // ���ݱ߳��򣬼Ӽ������ĵ�Ԫ������Ϊ�ӣ�����Ϊ��
    for (int i = 0; i < 4; i++) {
        int edgeIDi = edgeID[i];
        if (edgeIDi == -1) {
            continue;
        }
        for (int jVar = 0; jVar < 4; jVar++) {
            elementField_device.Flux[jVar][iElement] += edgeSign[i] * edgeField_device.Flux[jVar][edgeIDi];
        }
    }
}

// ������ճͨ��
void edge_convection_flux_device() {
    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementSoA& element_device = solver->element_device;
    GPU::ElementFieldSoA& elementField_device = solver->elementField_device;
    GPU::EdgeSoA& edge_device = solver->edge_device;
    GPU::EdgeFieldSoA& edgeField_device = solver->edgeField_device;
    GPU::BoundarySetMap& boundary_device = solver->boundary_device;

    int flux_conservation_scheme = GlobalPara::inviscid_flux_method::flux_conservation_scheme;
    int inviscid_flux_method_flag_reconstruct = GlobalPara::inviscid_flux_method::flag_reconstruct;
    myfloat* inf_ruvp = nullptr;
    myfloat* inlet_ruvp = nullptr;
    myfloat* outlet_ruvp = nullptr;
    cudaMalloc(&inf_ruvp, 4 * sizeof(myfloat));
    cudaMalloc(&inlet_ruvp, 4 * sizeof(myfloat));
    cudaMalloc(&outlet_ruvp, 4 * sizeof(myfloat));
    cudaMemcpy(inf_ruvp, GlobalPara::boundaryCondition::_2D::inf::ruvp, 4 * sizeof(myfloat), cudaMemcpyHostToDevice);
    cudaMemcpy(inlet_ruvp, GlobalPara::boundaryCondition::_2D::inlet::ruvp, 4 * sizeof(myfloat), cudaMemcpyHostToDevice);
    cudaMemcpy(outlet_ruvp, GlobalPara::boundaryCondition::_2D::outlet::ruvp, 4 * sizeof(myfloat), cudaMemcpyHostToDevice);
    myfloat gamma = GlobalPara::constant::gamma;
    myfloat rcpcv = GlobalPara::constant::R;
    // ����Ϊ˫��շ���ĳ�Ա����
    CBoundaryDoubleShockReflect* cbdsr = CBoundaryDoubleShockReflect::getInstance();
    myfloat dsr_shock_x = cbdsr->get_shock_x();
    myfloat dsr_shock_y = cbdsr->get_shock_y();
    myfloat dsr_shock_speed = cbdsr->get_shock_speed();
    myfloat dsr_t_RK = cbdsr->get_t_plus_dt();
    const myfloat _sqrt3 = sqrt(3.0);

    int block_size = GPU::MY_BLOCK_SIZE / 2;
    int grid_size = (edge_device.num_edge + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    edge_convection_flux_device_kernel <<<grid, block>>> (element_device, elementField_device, edge_device, edgeField_device, boundary_device,
        flux_conservation_scheme, inviscid_flux_method_flag_reconstruct, inf_ruvp, inlet_ruvp, outlet_ruvp, gamma, rcpcv, dsr_shock_x, dsr_shock_y, dsr_shock_speed,
        dsr_t_RK, _sqrt3);
    cudaFree(inf_ruvp);
    cudaFree(inlet_ruvp);
    cudaFree(outlet_ruvp);
    getLastCudaError("edge_convection_flux_device failed.");
}

void add_edge_flux_to_element_flux_device() {
    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementSoA& element_device = solver->element_device;
    GPU::ElementFieldSoA& elementField_device = solver->elementField_device;
    GPU::EdgeSoA& edge_device = solver->edge_device;
    GPU::EdgeFieldSoA& edgeField_device = solver->edgeField_device;
    int block_size = GPU::MY_BLOCK_SIZE;
    int grid_size = (element_device.num_element + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    add_edge_flux_to_element_flux_kernel <<<grid, block >>> (element_device, elementField_device, edge_device, edgeField_device);
    getLastCudaError("add_edge_flux_to_element_flux_device failed.");
}

void GPU::Space::Flux::calculateFluxDevice_2(ElementSoA& element_device, EdgeSoA& edge_device, ElementFieldSoA& elementField_device) {
    /*
    https://forums.developer.nvidia.com/t/synchronization-between-kernel-calls/23336
    ���κ˺�������֮�������cudaDeviceSynchronize().
    */

    clear_element_flux_device();
    edge_convection_flux_device();
    if (GlobalPara::physicsModel::equation == _EQ_NS) GPU::Space::edge_viscous_flux_device();
    add_edge_flux_to_element_flux_device();

    getLastCudaError("GPU::Space::Flux::calculateFluxDevice_2 failed.");
}
