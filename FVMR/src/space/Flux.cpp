#include "Flux.h"
#include "../legacy/FVM_2D.h"
#include "../math/Math.h"
#include "convection/Convection.h"
#include "../global/StringProcessor.h"
#include "restrict/Restrict.h"
#include "../output/LogWriter.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"
#include "../global/CExit.h"
#include "../solvers/SolverDataGetter.h"
#include "../solvers/GPUSolver2.h"
#include "viscous_flux/ViscousFlux.h"

/*
�ò��ִ������ڼ���ͨ����
��Ҫ�ñ���edgeField�洢���ߵ���ֵͨ�����ñ���������ⲿ�ļ���ȡ������������host��device֮�佻����
����֮������host���㣬ֻ������edgeField_host������device���㣬ֻ������edgeField_device
����host���룬Ϊ�˱����޸���;���к����������ȫ�ֱ����ṩ�����ⲿ�����Ľӿ�
����device���룬���ֱ�ӷ���ȫ�ֺ�����
*/
inline bool isUpStreamOfShock_atBoundary_static_inline(myfloat x, myfloat y, myfloat _shock_x, myfloat _shock_y, myfloat _shock_speed, myfloat _t_RK, const myfloat _sqrt3) {
    // CBoundaryDoubleShockReflect�ĳ�Ա�����ľ�̬�汾
    myfloat right = _shock_x + (y + _shock_speed * _t_RK) / _sqrt3;
    if (x < right) return true;
    else return false;
}

inline void get_U_reconstruct_singleValue(myfloat x, myfloat y, myfloat& U_dist, myint i_e, const myfloat* x_e, const myfloat* y_e, const myfloat* U_e, const myfloat* Ux_e, const myfloat* Uy_e, int flag_reconstruct) {
    if (flag_reconstruct == _REC_constant) {
        U_dist = U_e[i_e];
    }
    else if (flag_reconstruct == _REC_linear) {
        U_dist = U_e[i_e] + Ux_e[i_e] * (x - x_e[i_e]) + Uy_e[i_e] * (y - y_e[i_e]);
    }
    else {
        LogWriter::logAndPrintError("invalid flag_reconstruct. @get_U_reconstruct_singleValue.\n");
        U_dist = U_e[i_e];
    }
}

void get_Uvector_reconstruct_2(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, myfloat x, myfloat y, myint iElement, myfloat* U_dist, int flag_reconstruct, myfloat gamma) {
    // ���ݵ�Ԫ�ֲ��������õ�ĳ��Uֵ
    for (int i = 0; i < 4; i++) {
        get_U_reconstruct_singleValue(x, y, U_dist[i], iElement, element.xy[0], element.xy[1], elementField.U[i], elementField.Ux[i], elementField.Uy[i], flag_reconstruct);
    }

    // �������쳣�������ع�
    myfloat ruvp[4]{};
    U2NITS::Math::U2ruvp_host(U_dist, ruvp, gamma);
    if (U2NITS::Space::Restrict::outOfRange(ruvp)) {
        for (int i = 0; i < 4; i++) {
            U_dist[i] = elementField.U[i][iElement];
        }
    }
}

void get_UR_wallNonViscous(const myfloat* U_L, myfloat* U_R, myfloat nx, myfloat ny) {
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

void get_UR_wall_adiabat(const myfloat* U_L, myfloat* U_R, myfloat nx, myfloat ny) {
    // �޻��ƾ��ȱ��棬���������ٶȶ�ȡ��
    U_R[0] = U_L[0];
    U_R[1] = - U_R[1];
    U_R[2] = - U_R[2];
    U_R[3] = U_L[3];
}

void get_UR_farField(const myfloat* U_L, const myfloat* ruvp_inf, myfloat* U_R, myfloat nx, myfloat ny, myfloat gamma) {
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
    U2NITS::Math::U2ruvp_host(U_L, ruvp_e, gamma);
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
        // -as < qn < 0 ע��C++�в�������дС�ںţ�Ҫ��� -as<qn && qn<0
        // ���������
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
    // ���rho��p�Ƿ��쳣
    if (U2NITS::Space::Restrict::outOfRange(ruvp0)) {
        LogWriter::logAndPrintError("Out of range @get_UR_farField \n");
        CExit::saveAndExit(-1);
    }
    // ����ԭʼ����ruvp0�����غ���U_R
    U2NITS::Math::ruvp2U_host(ruvp0, U_R, gamma);
}

void get_UR_inner_and_periodic(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, GPU::EdgeSoA& edge, myint iEdge, myint iElementR, myfloat x_edge, myfloat y_edge, myfloat* U_R, int inviscid_flux_method_flag_reconstruct, myfloat gamma) {
    if (edge.setID[iEdge] != -1) {
        // ���ڱ߽� Ӧʹ�öԳƵ�edge���������Uֵ
        int ID = edge.ID[iEdge];
        int ID_pair = edge.periodicPair[ID];
        if (ID_pair < 0 || ID_pair >= edge.num_edge) {
            LogWriter::logAndPrintError("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel\n");
            exit(-1);
        }

        myfloat x_edge_pair = edge.xy[0][ID_pair];
        myfloat y_edge_pair = edge.xy[1][ID_pair];
        get_Uvector_reconstruct_2(element, elementField, x_edge_pair, y_edge_pair, iElementR, U_R, inviscid_flux_method_flag_reconstruct, gamma);
    }
    else {
        // �ڲ��߽�
        get_Uvector_reconstruct_2(element, elementField, x_edge, y_edge, iElementR, U_R, inviscid_flux_method_flag_reconstruct, gamma);
    }

}

void get_flux_RiemannSolve(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny, const myfloat length, myfloat* flux, const int flux_conservation_scheme, myfloat gamma, myfloat rcpcv) {
    // ���������������UL��UR���flux
    myfloat faceNormal[2]{ nx,ny };

    switch (flux_conservation_scheme) {// GlobalPara::inviscid_flux_method::flux_conservation_scheme
    case _SOL_LocalLaxFriedrichs:
        U2NITS::Space::LocalLaxFriedrichs(UL, UR, nx, ny, length, flux, gamma);
        break;
    case _SOL_Roe:
        U2NITS::Space::ConvectRoeCommon2d(UL, UR, faceNormal, length, flux, gamma, rcpcv);
        break;
    default:
        LogWriter::logAndPrintError("invalid scheme @get_flux_RiemannSolve");
        exit(-1);
        break;
    }
}

inline void fill_matrix(myfloat* mat[], int nrow, myint ncol, myfloat _Val) {
    /*
    ���룺matͨ��Ϊ4xn������������Ϊmyfloat* mat[4]��4��ָ�����ָ��һ��һά���飬��˸��з��������䣬Ҫ�ֱ�memset
    n�ϴ�(���ܳ���21��)�������myint

    memset���ֽ�Ϊ��λ���и�ֵ������ֵ�ķ�Χ��0x00��0xFF: memset(mat[i], 0, ncol * sizeof(myfloat));
    https://zhuanlan.zhihu.com/p/479101175
    int��ֵΪ0������_Val��0����
    int��ֵΪ1������_Val����ֱ�Ӹ�1����Ϊʵ���ϱ�ʾ������Ϊ0x01010101��û�бȽϺõķ���
    float��ֵΪ0��

    �����Ƽ���std::fill���ܸ�����ֵ����memset��ȫ����loop��
    https://iq.opengenus.org/set-array-to-0-in-cpp/
    */
    for (int i = 0; i < nrow; i++) {
        std::fill(mat[i], mat[i] + ncol, _Val);
    }
}

inline void clear_element_flux() {
    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementFieldSoA& elementField_host = solver->elementField_host;
    // ��Ԫ��ֵͨ��Flux����
    fill_matrix(elementField_host.Flux, 4, elementField_host.num, 0);
}

void edge_convection_flux_kernel(GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, GPU::EdgeFieldSoA& edgeField_host, GPU::BoundarySetMap& boundary_host,
    int flux_conservation_scheme, int inviscid_flux_method_flag_reconstruct, myfloat* inf_ruvp, myfloat* inlet_ruvp, myfloat* outlet_ruvp, myfloat gamma, myfloat rcpcv,
    myfloat dsr_shock_x, myfloat dsr_shock_y, myfloat dsr_shock_speed, myfloat dsr_t_RK, const myfloat _sqrt3, myint iEdge) {
    // �����GPUʱ������Ҫ������
    int setID = edge_host.setID[iEdge];// �ڲ�edge�������κ�set�����setIDΪ-1 setID�ĳ�ʼ����readContinueFile
    int bType = -1;
    if (setID != -1) {
        // ��iedge��edge�ı߽�����
        bType = boundary_host.type[setID - 1];
    }

    myfloat nx = edge_host.normal[0][iEdge];
    myfloat ny = edge_host.normal[1][iEdge];
    myfloat x_edge = edge_host.xy[0][iEdge];
    myfloat y_edge = edge_host.xy[1][iEdge];
    myfloat length = edge_host.length[iEdge];
    myint iElementL = edge_host.elementL[iEdge];
    myint iElementR = edge_host.elementR[iEdge];
    bool isUpStream_doubleShockReflect = isUpStreamOfShock_atBoundary_static_inline(x_edge, y_edge, dsr_shock_x, dsr_shock_y, dsr_shock_speed, dsr_t_RK, _sqrt3);

    myfloat U_L[4]{};
    myfloat U_R[4]{};
    myfloat flux[4]{};
    get_Uvector_reconstruct_2(element_host, elementField_host, x_edge, y_edge, iElementL, U_L, inviscid_flux_method_flag_reconstruct, gamma);// get U_L

    switch (bType) {
    case _BC_symmetry:
        // �Գơ���ŷ�����̣��൱����ճ�̱�
        get_UR_wallNonViscous(U_L, U_R, nx, ny);
        break;
    case _BC_wall_nonViscous:
        // ��ճ�̱�
        get_UR_wallNonViscous(U_L, U_R, nx, ny);
        break;
    case _BC_wall_adiabat:
        // �޻��ƾ���
        get_UR_wall_adiabat(U_L, U_R, nx, ny);
        break;
    case _BC_inlet:
        // ���
        get_UR_farField(U_L, inlet_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_outlet:
        // ����
        get_UR_farField(U_L, outlet_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_inf:
        // Զ��
        get_UR_farField(U_L, inf_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_doubleShockReflect:
        // ˫��շ���
        if (isUpStream_doubleShockReflect) {
            get_UR_farField(U_L, inlet_ruvp, U_R, nx, ny, gamma);
        }
        else {
            get_UR_farField(U_L, outlet_ruvp, U_R, nx, ny, gamma);
        }
        break;

    default:// �ڲ���bType=-1���߽磺bTypeȡ_BC_periodic_0��_BC_periodic_9����6100-6109
        if (iElementR != -1) {
            // ���ں��ڲ� ͳһ����
            get_UR_inner_and_periodic(element_host, elementField_host, edge_host, iEdge, iElementR, x_edge, y_edge, U_R, inviscid_flux_method_flag_reconstruct, gamma);
        }
    }

    get_flux_RiemannSolve(U_L, U_R, nx, ny, length, flux, flux_conservation_scheme, gamma, rcpcv);

    for (int j = 0; j < 4; j++) {
        edgeField_host.Flux[j][iEdge] = flux[j];
    }

}

// ������ճͨ�� created: 2024-05-01, tgl
void edge_convection_flux() {
    /*
    2023-05-05
    Ϊ�˸����GPU���룬��Ҫ����ǰ����
    ��ȫ�ֱ�����ɲ�����ʽ
    */

    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementSoA& element_host = solver->element_host;
    GPU::ElementFieldSoA& elementField_host = solver->elementField_host;
    GPU::EdgeSoA& edge_host = solver->edge_host;
    GPU::EdgeFieldSoA& edgeField_host = solver->edgeField_host;
    GPU::BoundarySetMap& boundary_host = solver->boundary_host;

    int flux_conservation_scheme = GlobalPara::inviscid_flux_method::flux_conservation_scheme;
    int inviscid_flux_method_flag_reconstruct = GlobalPara::inviscid_flux_method::flag_reconstruct;
    myfloat* inf_ruvp = GlobalPara::boundaryCondition::_2D::inf::ruvp;
    myfloat* inlet_ruvp = GlobalPara::boundaryCondition::_2D::inlet::ruvp;
    myfloat* outlet_ruvp = GlobalPara::boundaryCondition::_2D::outlet::ruvp;
    myfloat gamma = GlobalPara::constant::gamma;
    myfloat rcpcv = GlobalPara::constant::R;
    // ����Ϊ˫��շ���ĳ�Ա����
    CBoundaryDoubleShockReflect* cbdsr = CBoundaryDoubleShockReflect::getInstance();
    myfloat dsr_shock_x = cbdsr->get_shock_x();
    myfloat dsr_shock_y = cbdsr->get_shock_y();
    myfloat dsr_shock_speed = cbdsr->get_shock_speed();
    myfloat dsr_t_RK = cbdsr->get_t_plus_dt();
    const myfloat _sqrt3 = sqrt(3.0);

    // ÿ���߼�����ճͨ��
    for (myint iEdge = 0; iEdge < edge_host.num_edge; iEdge++) {
        edge_convection_flux_kernel(element_host, elementField_host, edge_host, edgeField_host, boundary_host,
            flux_conservation_scheme, inviscid_flux_method_flag_reconstruct, inf_ruvp, inlet_ruvp, outlet_ruvp, gamma, rcpcv, dsr_shock_x, dsr_shock_y, dsr_shock_speed,
            dsr_t_RK, _sqrt3, iEdge);
        
    }

}


inline void add_edge_flux_to_element_flux() {
    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementSoA& element_host = solver->element_host;
    GPU::ElementFieldSoA& elementField_host = solver->elementField_host;
    GPU::EdgeSoA& edge_host = solver->edge_host;
    GPU::EdgeFieldSoA& edgeField_host = solver->edgeField_host;

    // �Ӽ������ڵ�Ԫ
    for (myint iElement = 0; iElement < element_host.num_element; iElement++) {

        // �μ�Gradient
        // ��ǰ��Ԫ�����edgeID��edge�����������ⵥԪID
        myfloat volumeC = element_host.volume[iElement];// ���
        // ��ʼ��edgeID �������鷶Χ�ķǷ�ֵ��Ϊ-1
        int edgeID[4]{ -1,-1,-1,-1 };// ��iElement��element�ĵ�i���ߵ�ID
        for (int i = 0; i < 4; i++) {
            edgeID[i] = element_host.edges[i][iElement];
            if (!edge_host.has(edgeID[i])) {
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
            if (iElement != edge_host.elementR[edgeIDi]) {
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
                elementField_host.Flux[jVar][iElement] += edgeSign[i] * edgeField_host.Flux[jVar][edgeIDi];
            }
        }
    }
}

void U2NITS::Space::Flux::calculateFluxHost(GPU::ElementSoA& element_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {
    /*
    20240421: �ڱ����ɺ���calculate_flux_1��ͬʱ�����º���calculate_flux_2��inline�����˺������õĿ���
    �������ܣ����㵥Ԫ��ֵͨ�������ȱ������ߣ�����������غ�����Ȼ�������������õ���ͨ�������ֱ�Ӽ������ҵ�Ԫ
    ��֮ǰ����������
    ����õ���Flux��ֱ�ӼӼ����ھӵ�Ԫ�������ȴ���edgeField_host
    �ŵ�1������ѭ�����һ��edgeѭ����һ��elementѭ���������ڲ���
    �ŵ�2����Flux����edgeField_host�����ڲ���ճ��ͨ����ⲿ�֣�Ϊ�պ��ճ����׼��

    �ȱ����ߣ�������ճͨ�������edgeField_host
    Ȼ������ߣ�����ճ��ͨ������edgeField_host��flux���޸�
    Ȼ�������Ԫ��ȡedgeField_host��ֵ���Ӽ�������
    */
    clear_element_flux();// element��ֵͨ������
    edge_convection_flux();// edge��ճͨ��
    if (GlobalPara::physicsModel::equation == _EQ_NS) U2NITS::Space::edge_viscous_flux();// edgeճ��ͨ��
    add_edge_flux_to_element_flux();// edge��ֵͨ���Ӽ���element��ֵͨ��

}



