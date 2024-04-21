#include "Flux.h"
#include "../FVM_2D.h"
#include "../math/Math.h"
#include "convection/Convection.h"
#include "../global/StringProcessor.h"
#include "restrict/Restrict.h"
#include "../output/LogWriter.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"


void U2NITS::Space::Flux::calculateFluxHost(ElementSoA& element_host, EdgeSoA& edge_host, FieldSoA& elementField_host /*, std::map<int, int>& edge_periodic_pair*/ ) {
    /*
    �ο���void Solver_2D::calFlux()
    TODO:
    �������ڱ߽�Ĵ�����Ҫ�޸ġ����֧��10�����ڱ߽磬���и��࣬�ᰴinner���㣬���´���
    */

    // ��ʼ����Ԫ��ֵͨ�����ֲ�����
    for (int iElement = 0; iElement < element_host.num_element; iElement++) {
        for (int jValue = 0; jValue < 4; jValue++) {
            elementField_host.Flux[jValue][iElement] = 0;//��Ԫ��ֵͨ��Flux���㣬Ϊ����Ӽ���׼��
            // f->elements[ie].deltaeig = 0;//ÿһ��deltaeig���� ���Roe��ʽ Ŀǰ����Ҫ
            // �ݶ��Ѿ����㣬��˲���Ҫ�ٴμ���
        }
    }

    // ÿ���߼�����ճͨ����Ȼ����ݷ���ֱ�Ӽ������൥Ԫ��Flux�����б߱��������е�Ԫ��FluxҲ�ͼ��������
    for (int iedge = 0; iedge < edge_host.num_edge; iedge++) {

        FVM_2D* f = FVM_2D::getInstance();
        int elementL = edge_host.elementL[iedge];
        int elementR = edge_host.elementR[iedge];
        int setID = edge_host.setID[iedge];// �ڲ�edge�������κ�set�����setIDΪ-1 setID�ĳ�ʼ����readContinueFile
        int boundaryNum = f->boundaryManager.boundaries.size();
        int bType = -1;
        if (setID != -1) {
            // ��iedge��edge�ı߽����� = boundaries[setID - 1].type������setID = edge_host.setID[iedge]
            bType = f->boundaryManager.boundaries[setID - 1].type;
        }
        double flux[4]{};


        //if (iedge == 562) {
        //    //std::cout << "ElementL=" << elementL << ", x_edge=" << x_edge << ", y_edge=" << y_edge << ", ";
        //    std::cout << "edge=" << iedge << ", ";
        //    //std::cout << "ruvp_L=" << StringProcessor::doubleArray_2_string(ruvp_L, 4) << std::endl;
        //    std::cout << "(U2NITS::Space::Flux::getEdgeFlux_farField)\n";
        //}

        switch (bType) {
        case _BC_symmetry:
            // �Գơ���ŷ�����̣��൱����ճ�̱�
            getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iedge, flux);
            break;
        case _BC_wall_nonViscous:
            // ��ճ�̱�
            getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iedge, flux);
            break;
        case _BC_wall_adiabat:
            // �޻��ƾ���
            getEdgeFlux_wall_adiabat(element_host, elementField_host, edge_host, iedge, flux);
            break;
        case _BC_inlet:
            // ���
            getEdgeFlux_farField_3(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
            break;
        case _BC_outlet:
            // ����
            getEdgeFlux_farField_3(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
            break;
        case _BC_inf:
            // Զ��
            getEdgeFlux_farField_3(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::inf::ruvp, flux);
            break;
        case _BC_doubleShockReflect:
            // ˫��շ���
            getEdgeFlux_farField_doubleShockReflect(element_host, elementField_host, edge_host, iedge, flux);

            break;
        default:// �ڲ���bType=-1���߽磺bTypeȡ_BC_periodic_0��_BC_periodic_9����6100-6109
            if (elementR != -1) {
                // ���ں��ڲ� ͳһ����
                getEdgeFlux_inner_and_periodic(element_host, elementField_host, edge_host, iedge, flux);
            }
        }

        // ���µ�Ԫ��ֵͨ��
        for (int j = 0; j < 4; j++) {
            elementField_host.Flux[j][elementL] += flux[j];
            // �ڲ��߽�
            if (bType == -1) {
                // �˴�����elementR==-1����ֹ���ڱ߽���2��
                elementField_host.Flux[j][elementR] -= flux[j];
            }
        }
    }// end for iedge
}

REAL U2NITS::Space::Flux::calculateLambdaFlux(REAL edgeU[4], REAL edgeN[2], REAL gamma, REAL length) {

    REAL LambdaC = 0;
    for (int ie = 0; ie < 3; ie++) {
        //eU
        auto u = edgeU[1] / edgeU[0];
        auto v = edgeU[2] / edgeU[0];
        auto E = edgeU[3] / edgeU[0];
        auto eabs = abs(u * edgeN[0] + v * edgeN[1]);//|euv��en|
        //ec
        auto V2 = u * u + v * v;
        REAL& rho = edgeU[0];
        REAL p = 0.5 * rho * (gamma - 1) * (2 * E - V2);
        REAL ec = sqrt(gamma * p / rho);
        REAL& dl = length;
        LambdaC += (eabs + ec) * dl;
    }
    return LambdaC;


}

void U2NITS::Space::Flux::getEdgeFlux_wallNonViscous(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, long idx, double* flux) {
    // ������idx�൱��iedge
    // �ú����ֱ����ߵ���������U�ļ��ޣ�Ȼ�����flux

    REAL nx = edge.normal[0][idx];// normal�ѳ�ʼ��
    REAL ny = edge.normal[1][idx];
    REAL x_edge = edge.xy[0][idx];
    REAL y_edge = edge.xy[1][idx];
    REAL length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // ����U_L����edge U������
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementL, U_L);

    // ����U_R �öԳ���
    REAL U_R[4]{};
    REAL uL = U_L[1] / U_L[0];
    REAL vL = U_L[2] / U_L[0];
    REAL uR = uL - 2 * (uL * nx + vL * ny) * nx;
    REAL vR = vL - 2 * (uL * nx + vL * ny) * ny;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    // ����flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getEdgeFlux_wall_adiabat(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, long idx, double* flux) {
    // �̱ڱ߽硣������idx�൱��iedge

    REAL nx = edge.normal[0][idx];// normal�ѳ�ʼ��
    REAL ny = edge.normal[1][idx];
    REAL x_edge = edge.xy[0][idx];
    REAL y_edge = edge.xy[1][idx];
    REAL length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // ����U_L����edge U������
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementL, U_L);

    // ����U_R �öԳ���
    REAL U_R[4]{};
    REAL uL = U_L[1] / U_L[0];
    REAL vL = U_L[2] / U_L[0];
    // ��������ϵ�� UnL = uL * nx + vL * ny��VnL = -uL * ny + vL * nx
    // UnR = - UnL��VnR = -VnL
    // ԭ����ϵ�� uR = unR*nx-vnR*ny�� vR=unR*ny+vnR*nx
    // ��� uR = -uL ������ֱ��ȡ��

    // ��UNITS��BC_WALL.f90���㷨���ٶ�ȡ�������ǽ����Ȼ����
    // ���� �����ͼ��û������ճ��
    //REAL uR = uL - 2 * (uL * nx + vL * ny) * nx;
    //REAL vR = vL - 2 * (uL * nx + vL * ny) * ny;
    REAL uR = -uL;
    REAL vR = -vL;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    // ����flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getEdgeFlux_farField(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, long idx, REAL* ruvp_inf, REAL* flux) {
    
    constexpr bool useNew = true;
    if (useNew) {
        getEdgeFlux_farField_3(element, elementField, edge, idx, ruvp_inf, flux);
        return;
    }

    
    // �����Ǿɰ����
    REAL nx = edge.normal[0][idx];// normal�ѳ�ʼ��
    REAL ny = edge.normal[1][idx];
    REAL x_edge = edge.xy[0][idx];
    REAL y_edge = edge.xy[1][idx];
    REAL length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // ����U_L����edge U������
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementL, U_L);
    // ����ruvp_L
    REAL ruvp_L[4]{};
    U2NITS::Math::U2ruvp_host(U_L, ruvp_L, GlobalPara::constant::gamma);
    // ����Ƿ�ɢ
    if (Restrict::outOfRange(ruvp_L)) {
        std::stringstream ss;
        ss << "ElementL=" << iElementL << ", x_edge=" << x_edge << ", y_edge=" << y_edge << ", ";
        ss << "edge=" << idx << ", ";
        ss << "ruvp_L=" << StringProcessor::doubleArray_2_string(ruvp_L, 4) << std::endl;
        ss << "(U2NITS::Space::Flux::getEdgeFlux_farField_kernel)\n";
        LogWriter::logAndPrintError(ss.str());
    }

    //����Զ���߽���������ruvp_L
    Space::Flux::modify_ruvpL_farField(nx, ny, ruvp_L, ruvp_inf);
    U2NITS::Math::ruvp2U_host(ruvp_L, U_L, GlobalPara::constant::gamma);
    // ����flux
    RiemannSolve(U_L, U_L, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getEdgeFlux_farField_2(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long iEdge, REAL* ruvp_inf, REAL* flux) {
    /*
    20240412 ԭ����Զ���߽�������
    20240415 �о������⣬������Զ����Ϊ�оݡ���μ�UNITs: BC/BC_far.f90: subroutine BC_farbound
    �½�����farField_3�����޸�
    */

    if (!edge_host.has(iEdge)) {
        LogWriter::logAndPrintError("iEdge" + std::to_string(iEdge) + " out of range.\n");
        exit(-1);
    }
    // Զ������ inf
    real nx = edge_host.normal[0][iEdge];
    real ny = edge_host.normal[1][iEdge];
    real gamma = GlobalPara::constant::gamma;
    real rho_inf = ruvp_inf[0];
    real u_inf = ruvp_inf[1];
    real v_inf = ruvp_inf[2];
    real p_inf = ruvp_inf[3];
    real un_inf = u_inf * nx + v_inf * ny;
    real a2_inf = gamma * p_inf / rho_inf;
    real Ma2_inf = (u_inf * u_inf + v_inf * v_inf) / a2_inf;
    real Man2_inf = (un_inf * un_inf) / a2_inf;

    // �ڵ�Ԫ���� domain �����ű�d��ʾdomain
    real U_d[4]{};
    integer element_inner = edge_host.elementL[iEdge];
    for (int i = 0; i < 4; i++) {
        U_d[i] = elementField_host.U[i][element_inner];
    }
    real ruvp_d[4]{};
    Math::U2ruvp_host(U_d, ruvp_d, gamma);
    real rho_d = ruvp_d[0];
    real u_d = ruvp_d[1];
    real v_d = ruvp_d[2];
    real p_d = ruvp_d[3];
    real un_d = u_d * nx + v_d * ny;
    real a2_d = gamma * p_d / rho_d;
    real Ma2_d = (u_d * u_d + v_d * v_d) / a2_d;
    real Man2_d = (un_d * un_d) / a2_d;

    // �߽���� boundary ����
    bool isFarfieldSuperSonic = (Ma2_inf > 1);// �г����²�����Ϊ�˲����ĸ����ʺ���Ϊ�о�
    bool isFarfieldNormalSuperSonic = (Man2_inf > 1);
    bool isDomainSuperSonic = (Ma2_d > 1);
    bool isDomainNormalSuperSonic = (Man2_d > 1);
    bool isFarfieldInlet = (un_inf < 0);
    bool isDomainInlet = (un_d < 0);
    bool isSuperSonic = isDomainNormalSuperSonic;
    bool isInlet = isDomainInlet;
    real U_b[4]{};// U at boundary edge
    if (isSuperSonic) {
        if (isInlet) {
            Math::ruvp2U_host(ruvp_inf, U_b, gamma);
        }
        else {
            for (int i = 0; i < 4; i++) {
                U_b[i] = U_d[i];
            }
        }
    }
    else {
        // �ο� CFDPA P264 reference state
        real rho_a = rho_inf;// �ű��Ϊa�������б���һ�£�����У��
        real u_a = u_inf;
        real v_a = v_inf;
        real p_a = p_inf;

        real rho0 = rho_d;
        real c02 = gamma * p_d / rho_d;
        real c0 = sqrt(c02);
        real rho0c0 = rho0 * c0;

        if (isInlet) {
            real p_b = 0.5 * (p_a + p_d - rho0c0 *
                (nx*(u_a-u_d)+ny*(v_a-v_d)) );
            real pa_pb = p_a - p_b;
            real rho_b = rho_a + (-pa_pb) / c02;
            real u_b = u_a - nx * pa_pb / rho0c0;
            real v_b = v_a - ny * pa_pb / rho0c0;
            real ruvp_b[4]{ rho_b,u_b,v_b,p_b };
            Math::ruvp2U_host(ruvp_b, U_b, gamma);
        }
        else {
            real p_b = p_a;
            real pd_pb = p_d - p_b;
            real rho_b = rho_d + (-pd_pb) / c02;
            real u_b = u_d + nx * pd_pb / rho0c0;
            real v_b = v_d + ny * (pd_pb) / rho0c0;
            real ruvp_b[4]{ rho_b,u_b,v_b,p_b };
            Math::ruvp2U_host(ruvp_b, U_b, gamma);
        }
    }
    // ����flux
    real length = edge_host.length[iEdge];
    RiemannSolve(U_d, U_b, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getEdgeFlux_farField_3(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, integer iEdge, REAL* ruvp_inf, REAL* flux) {
    /*
    ����UNITs���������μ�
    �����������Ϊ���������ڵ����Ʊ߽��Ĺ���ֵ��
    �ٽ����߽����ֵ�������ֵ���������ݹ�ϵ
    ����߽����٣������õ�ѹ����ֵ��
    */

    if (!edge_host.has(iEdge)) {
        LogWriter::logAndPrintError("iEdge" + std::to_string(iEdge) + " out of range.\n");
        exit(-1);
    }
    // Զ��ֵ
    real nx = edge_host.normal[0][iEdge];
    real ny = edge_host.normal[1][iEdge];
    real gamma = GlobalPara::constant::gamma;
    real ga1 = gamma - 1;
    real rho_f = ruvp_inf[0];
    real u_f = ruvp_inf[1];
    real v_f = ruvp_inf[2];
    real p_f = ruvp_inf[3];
    
    // �ڵ�ֵ����߽��������ֵ �˴��õ�Ԫ����ֵ��Ϊ����ֵ
    real U_e[4]{};
    integer element_inner = edge_host.elementL[iEdge];
    for (int i = 0; i < 4; i++) {
        U_e[i] = elementField_host.U[i][element_inner];
    }
    real ruvp_e[4]{};
    Math::U2ruvp_host(U_e, ruvp_e, gamma);
    real rho_e = ruvp_e[0];
    real u_e = ruvp_e[1];
    real v_e = ruvp_e[2];
    real p_e = ruvp_e[3];

    // �������ݹ�ϵ
    real a2_f = gamma * p_f / rho_f;
    real af = sqrt(a2_f);// Զ������
    real Ma2_f = (u_f * u_f + v_f * v_f) / a2_f;// Զ�������ƽ��
    real ae = sqrt(gamma * p_e / rho_e);// �ڵ�����
    real qnf = u_f * nx + v_f * ny;// Զ�������ٶ�
    real qne = u_e * nx + v_e * ny;// �ڵ㷨���ٶ�
    real rf = qnf - 2.0 * af / ga1;// ��һ������ R2 = u0n - 2*a0n/(gamma-1)
    real re = qne + 2.0 * ae / ga1;// �ڶ������� u0n - 2*a0n/(gamma-1)
    real qn = 0.5 * (re + rf);
    real as = ga1 * (re - rf) / 4.0;
    real dnt = 0;// �����˶��ٶ�

    // �о� q<=-as����������ڣ�-as<q<0����������ڣ�0<=q<as�������ٳ��ڣ�as<=q�������ٳ���
    real rho0 = 1, u0 = 0, v0 = 0, p0 = 1;
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
        real son = p_f / pow(rho_f, gamma);// ��
        real qtx = u_f - qnf * nx;
        real qty = v_f - qnf * ny;

        rho0 = pow((as * as / son / gamma), 1.0 / ga1);// ����S=p/rho^gamma��rho=(p/S)^(1/gamma)=(a2/gamma/S)^(1/gamma)
        u0 = qtx + qn * nx;
        v0 = qty + qn * ny;
        p0 = as * as * rho0 / gamma;// p=rho*R*t=rho/gamma*gamma*R*t=rho/gamma*a2
    }
    else if (qn < as) {
        // 0 <= qn < as
        // �����ٳ��� 
        real son = p_e / pow(rho_e, gamma);
        real qtx = u_e - qne * nx;
        real qty = v_e - qne * ny;

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
    real ruvp0[4]{ rho0,u0,v0,p0 };
    // ����flux
    real U0[4]{};
    Math::ruvp2U_host(ruvp0, U0, gamma);
    if (Restrict::outOfRange(ruvp0)) {
        LogWriter::logAndPrintError("Out of range @U2NITS::Space::Flux::getEdgeFlux_farField_3 \n");
    }

    real length = edge_host.length[iEdge];
    RiemannSolve(U_e, U0, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);

}

void U2NITS::Space::Flux::modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, const REAL* ruvp_inf) {
    // �ú�������Զ���߽���������ruvp_L

    //�߽�Ԫ��Element_RΪnull�����nxnyһ������
    const REAL rho = ruvp[0];
    const REAL u = ruvp[1];
    const REAL v = ruvp[2];
    const REAL p = ruvp[3];
    const REAL rho_inf = ruvp_inf[0];
    const REAL u_inf = ruvp_inf[1];
    const REAL v_inf = ruvp_inf[2];
    const REAL p_inf = ruvp_inf[3];
    const REAL gamma = GlobalPara::constant::gamma;

    //����任
    const REAL rho_n = rho;
    const REAL u_n = u * nx + v * ny;//u cost + v sint
    const REAL v_n = -u * ny + v * nx;//�˴����޸ġ���������Ϊ�� - u sint + v cost
    const REAL p_n = p;
    const REAL rho_n_inf = rho_inf;
    const REAL p_n_inf = p_inf;
    const REAL u_n_inf = u_inf * nx + v_inf * ny;
    const REAL v_n_inf = -u_inf * ny + v_inf * nx;//�˴����޸�

    
    REAL a_n2 = gamma * p_n / rho_n;
    if (a_n2 < 0) {
        std::cout << "Error: a_n2 < 0\n";  
        
    }
    const REAL a_n = sqrt(a_n2);
    REAL a_n_inf2 = gamma * p_n_inf / rho_n_inf;
    if (a_n_inf2 < 0) {
        std::cout << "Error: a_n_inf2 < 0\n";  
        
    }
    const REAL a_n_inf = sqrt(a_n_inf2);

    REAL v_n_tmp;
    REAL u_n_tmp;
    REAL rho_n_tmp;
    REAL p_n_tmp;

    if (u_n_inf * u_n_inf < a_n_inf * a_n_inf) {//|Vn_inf|<|a_inf|��������
        REAL tmp_a = (gamma - 1.0) / 4.0 * (u_n - u_n_inf) + 0.5 * (a_n + a_n_inf);// 1/2 *( (ga-1)/2 * (u-uinf) + (a+ainf) )
        REAL tmp_a2 = tmp_a * tmp_a;
        if (u_n_inf <= 0.0) {//��������ڣ���3���߽�����
            REAL tmp_s_inf = p_n_inf / pow(rho_n_inf, gamma);//���ݱ߽�rho��p�����
            rho_n_tmp = pow(tmp_a2 / tmp_s_inf / gamma, 1.0 / (gamma - 1.0));//�߽�
            u_n_tmp = 0.5 * (u_n + u_n_inf + 2.0 * a_n / (gamma - 1.0) - 2.0 * a_n_inf / (gamma - 1.0));
            v_n_tmp = v_n_inf;//�߽�
            p_n_tmp = tmp_a2 * rho_n_tmp / gamma;//�߽�
        }
        else {//�����ٳ��ڣ���1���߽����� u_n_tmp
            REAL tmp_s = p_n / pow(rho_n, gamma);//�� //pow(x,y):��x<0��y������������x=0��y<=0�������ֽ������
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

void U2NITS::Space::Flux::getEdgeFlux_inner_and_periodic(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, long idx, REAL* flux) {


    //������е����ꡢ����任����������任���󡢱߷��߷���
    REAL nx = edge.normal[0][idx];// normal�ѳ�ʼ��
    REAL ny = edge.normal[1][idx];
    REAL x_edge = edge.xy[0][idx];
    REAL y_edge = edge.xy[1][idx];
    REAL length = edge.length[idx];
    long iElementL = edge.elementL[idx];
    long iElementR = edge.elementR[idx];

    // ����edge���괦��Uֵ(�ֱ������൥Ԫ�ķֲ���������)����Ϊ���ڱ߽���ڲ��߽��������
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementL, U_L);
    REAL U_R[4]{};
    if (edge.setID[idx] != -1) {
        // ���ڱ߽� Ӧʹ�öԳƵ�edge���������Uֵ
        int ID = edge.ID[idx];

        //// ��һ ��map
        //auto cIter = edge_periodic_pair.find(ID);
        //if (cIter == edge_periodic_pair.end()) {
        //    LogWriter::logAndPrintError("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic\n");
        //    exit(-1);
        //}
        //int ID_pair = cIter->second;

        // ���� ��edge�洢
        int ID_pair = edge.periodicPair[ID];
        if (ID_pair < 0 || ID_pair >= edge.num_edge) {
            LogWriter::logAndPrintError("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel\n");
            exit(-1);
        }

        REAL x_edge_pair = edge.xy[0][ID_pair];
        REAL y_edge_pair = edge.xy[1][ID_pair];
        Flux::getUByXYandElementID(element, elementField, x_edge_pair, y_edge_pair, iElementR, U_R);
    }
    else {
        // �ڲ��߽�
        Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementR, U_R);
    }


    // ����flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getUByXYandElementID(ElementSoA& element, FieldSoA& elementField, REAL x, REAL y, int elementID, REAL* U_dist) {
    // �ú������ݵ�Ԫ���ع�����������ĳ���괦��Uֵ
    REAL x_elementL = element.xy[0][elementID];
    REAL y_elementL = element.xy[1][elementID];
    // �����ع�
    if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_constant) {
        for (int jValue = 0; jValue < 4; jValue++) {
            U_dist[jValue] = elementField.U[jValue][elementID];
        }
    }
    // �����ع� �����ݶ�Ux��Uy�����(xpoint,ypoint)��_Uֵ
    else if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_linear) {
        for (int jValue = 0; jValue < 4; jValue++) {
            REAL& U_elementL = elementField.U[jValue][elementID];
            REAL& Ux_elementL = elementField.Ux[jValue][elementID];
            REAL& Uy_elementL = elementField.Uy[jValue][elementID];
            U_dist[jValue] = U_elementL + Ux_elementL * (x - x_elementL) + Uy_elementL * (y - y_elementL);
        }
        // �������쳣�������ع�
        real ruvp[4]{};
        Math::U2ruvp_host(U_dist, ruvp, GlobalPara::constant::gamma);
        if (Restrict::outOfRange(ruvp)) {
            for (int jValue = 0; jValue < 4; jValue++) {
                U_dist[jValue] = elementField.U[jValue][elementID];
            }
        }
    }
    else {
        LogWriter::logAndPrintError("implemented. @getUByXYandElementID_kernel.\n");
    }
}

void U2NITS::Space::Flux::getEdgeFlux_farField_doubleShockReflect(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, integer iEdge, REAL* flux) {
    /*
    ����߽�Ԫ��λ�����λ������Σ����ֱ�����ںͳ��ڵĲ�����ΪԶ���߽����
    */
    real edgex = edge_host.xy[0][iEdge];
    real edgey = edge_host.xy[1][iEdge];
    CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
    bool isUpStream = pCDS->isUpStreamOfShock_atBoundary(edgex, edgey, pCDS->get_t_plus_dt());
    if (isUpStream) {
        Space::Flux::getEdgeFlux_farField_3(element_host, elementField_host, edge_host, iEdge,
            GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
    }
    else {
        Space::Flux::getEdgeFlux_farField_3(element_host, elementField_host, edge_host, iEdge,
            GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
    }
}

void U2NITS::Space::Flux::RiemannSolve(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny, const REAL length, REAL* flux, const int scheme) {
    real faceNormal[2]{ nx,ny };
    real gamma = GlobalPara::constant::gamma;
    switch (scheme) {// GlobalPara::inviscid_flux_method::flux_conservation_scheme
    case _SOL_LocalLaxFriedrichs:
        U2NITS::Space::LocalLaxFriedrichs(UL, UR, nx, ny, length, flux, gamma);
        break;
    case _SOL_Roe:
        U2NITS::Space::ConvectRoeCommon2d(UL, UR, faceNormal, length, flux, gamma);
        break;
    default:
        break;
    }
}



