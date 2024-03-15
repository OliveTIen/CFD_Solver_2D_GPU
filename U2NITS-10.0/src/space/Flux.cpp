#include "Flux.h"
#include "../FVM_2D.h"
#include "../math/PhysicalConvertKernel.h"
#include "convection/Convection.h"



void U2NITS::Space::Flux::calculateFluxHost(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
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

        switch (bType) {
            //�Գơ���ŷ�����̣��൱����ճ�̱�
        case _BC_symmetry:
            Space::Flux::getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iedge, flux);
            break;
            //��ճ�̱�
        case _BC_wall_nonViscous:
            Space::Flux::getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iedge, flux);
            break;
            //�޻��ƾ���
        case _BC_wall_adiabat:
            // TODO: ʵ���޻��ƾ��ȱ߽������ļ���
            Space::Flux::getEdgeFlux_wall_adiabat(element_host, elementField_host, edge_host, iedge, flux);
            break;
            //���
        case _BC_inlet:
            Space::Flux::getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
            break;
            //����
        case _BC_outlet:
            Space::Flux::getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
            break;
            //Զ��
        case _BC_inf:
            Space::Flux::getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::inf::ruvp, flux);
            break;
        default:// �ڲ���bType=-1
            if (elementR != -1) {
                // ���ں��ڲ� ͳһ����
                Space::Flux::getEdgeFlux_inner_and_periodic(element_host, elementField_host, edge_host, edge_periodic_pair, iedge, flux);
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

void U2NITS::Space::Flux::getEdgeFlux_wallNonViscous(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux) {
    // ������idx�൱��iedge
    // �ú����ֱ����ߵ���������U�ļ��ޣ�Ȼ�����flux

    REAL nx = edge_host.normal[0][idx];// normal�ѳ�ʼ��
    REAL ny = edge_host.normal[1][idx];
    REAL x_edge = edge_host.xy[0][idx];
    REAL y_edge = edge_host.xy[1][idx];
    REAL length = edge_host.length[idx];
    long iElementL = edge_host.elementL[idx];

    // ����U_L����edge U������
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element_host, elementField_host, x_edge, y_edge, iElementL, U_L);

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

void U2NITS::Space::Flux::getEdgeFlux_wall_adiabat(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux) {
    // �̱ڱ߽硣������idx�൱��iedge

    REAL nx = edge_host.normal[0][idx];// normal�ѳ�ʼ��
    REAL ny = edge_host.normal[1][idx];
    REAL x_edge = edge_host.xy[0][idx];
    REAL y_edge = edge_host.xy[1][idx];
    REAL length = edge_host.length[idx];
    long iElementL = edge_host.elementL[idx];

    // ����U_L����edge U������
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element_host, elementField_host, x_edge, y_edge, iElementL, U_L);

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

void U2NITS::Space::Flux::getEdgeFlux_farField(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* ruvp_inf, REAL* flux) {
    // �������ܣ�����Զ���߽���������߽���ֵͨ�� �����

    REAL nx = edge_host.normal[0][idx];// normal�ѳ�ʼ��
    REAL ny = edge_host.normal[1][idx];
    REAL x_edge = edge_host.xy[0][idx];
    REAL y_edge = edge_host.xy[1][idx];
    REAL length = edge_host.length[idx];
    long iElementL = edge_host.elementL[idx];

    // ����U_L����edge U������
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element_host, elementField_host, x_edge, y_edge, iElementL, U_L);
    // ����ruvp_L
    REAL ruvp_L[4]{};
    U2NITS::Math::U2ruvp_host(U_L, ruvp_L, GlobalPara::constant::gamma);
    // ����Ƿ�ɢ
    bool is_divergent = false;
    if (ruvp_L[3] < 0) {
        is_divergent = true;
        std::cout << "p < 0, ";
    }
    if (ruvp_L[0] < 0) {
        is_divergent = true;
        std::cout << "rho < 0, ";
    }
    if (ruvp_L[0] > 900) {
        is_divergent = true;
        std::cout << "rho > 900, ";
    }
    if (is_divergent) {
        std::cout << "ElementLID=" << iElementL << ", x_edge=" << x_edge << ", y_edge=" << y_edge << ", ";
        std::cout << "(GPU::Flux::getEdgeFlux_farField)\n";
    }
    //����Զ���߽���������ruvp_L
    Space::Flux::modify_ruvpL_farField(nx, ny, ruvp_L, ruvp_inf);
    U2NITS::Math::ruvp2U_host(ruvp_L, U_L, GlobalPara::constant::gamma);
    // ����flux
    RiemannSolve(U_L, U_L, nx, ny, length, flux,
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

    //if (rho_n < 0)system("pause");
    REAL a_n2 = gamma * p_n / rho_n;
    if (a_n2 < 0) { std::cout << "Error: a_n2 < 0\n";  system("pause"); }
    const REAL a_n = sqrt(a_n2);
    REAL a_n_inf2 = gamma * p_n_inf / rho_n_inf;
    if (a_n_inf2 < 0) { std::cout << "Error: a_n_inf2 < 0\n";  system("pause"); }
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

void U2NITS::Space::Flux::getEdgeFlux_inner_and_periodic(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair, long idx, REAL* flux) {


    //������е����ꡢ����任����������任���󡢱߷��߷���
    REAL nx = edge_host.normal[0][idx];// normal�ѳ�ʼ��
    REAL ny = edge_host.normal[1][idx];
    REAL x_edge = edge_host.xy[0][idx];
    REAL y_edge = edge_host.xy[1][idx];
    REAL length = edge_host.length[idx];
    long iElementL = edge_host.elementL[idx];
    long iElementR = edge_host.elementR[idx];

    // ����edge���괦��Uֵ(�ֱ������൥Ԫ�ķֲ���������)
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element_host, elementField_host, x_edge, y_edge, iElementL, U_L);
    REAL U_R[4]{};

    if (edge_host.setID[idx] != -1) {
        // ���ڱ߽� Ӧʹ�öԳƵ�edge���������Uֵ
        //throw "Warning: periodic not implemented yet!\n";
        int ID = edge_host.ID[idx];
        auto cIter = edge_periodic_pair.find(ID);
        if (cIter == edge_periodic_pair.end()) {
            throw "Error: pair not found!\n";
        }
        else {
            int ID_pair = cIter->second;
            REAL x_edge_pair = edge_host.xy[0][ID_pair];
            REAL y_edge_pair = edge_host.xy[1][ID_pair];
            Flux::getUByXYandElementID(element_host, elementField_host, x_edge_pair, y_edge_pair, iElementR, U_R);
        }
    }
    else {
        // �ڲ��߽�
        Flux::getUByXYandElementID(element_host, elementField_host, x_edge, y_edge, iElementR, U_R);
    }

    // ����flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getUByXYandElementID(ElementSoA& element_host, FieldSoA& elementField_host, REAL x, REAL y, int elementID, REAL* U_dist) {
    // �ú������ݵ�Ԫ���ع�����������ĳ���괦��Uֵ
    REAL x_elementL = element_host.xy[0][elementID];
    REAL y_elementL = element_host.xy[1][elementID];
    // �����ع�
    if (GlobalPara::space::flag_reconstruct == _REC_constant) {
        for (int jValue = 0; jValue < 4; jValue++) {
            U_dist[jValue] = elementField_host.U[jValue][elementID];
        }
    }
    // �����ع� �����ݶ�Ux��Uy�����(xpoint,ypoint)��_Uֵ
    else if (GlobalPara::space::flag_reconstruct == _REC_linear) {
        for (int jValue = 0; jValue < 4; jValue++) {
            REAL& U_elementL = elementField_host.U[jValue][elementID];
            REAL& Ux_elementL = elementField_host.Ux[jValue][elementID];
            REAL& Uy_elementL = elementField_host.Uy[jValue][elementID];
            U_dist[jValue] = U_elementL + Ux_elementL * (x - x_elementL) + Uy_elementL * (y - y_elementL);
        }
    }
}

void U2NITS::Space::Flux::RiemannSolve(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny, const REAL length, REAL* flux, const int scheme) {
    real faceNormal[2]{ nx,ny };
    real gamma = GlobalPara::constant::gamma;
    switch (scheme) {
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



