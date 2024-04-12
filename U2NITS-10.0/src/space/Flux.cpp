#include "Flux.h"
#include "../FVM_2D.h"
#include "../math/Math.h"
#include "convection/Convection.h"
#include "../global/StringProcessor.h"
#include "restrict/Restrict.h"
#include "../output/LogWriter.h"


void U2NITS::Space::Flux::calculateFluxHost(ElementSoA& element_host, EdgeSoA& edge_host, FieldSoA& elementField_host /*, std::map<int, int>& edge_periodic_pair*/ ) {
    /*
    参考：void Solver_2D::calFlux()
    TODO:
    计算周期边界的代码需要修改。最多支持10对周期边界，若有更多，会按inner计算，导致错误
    */

    // 初始化单元数值通量、分布函数
    for (int iElement = 0; iElement < element_host.num_element; iElement++) {
        for (int jValue = 0; jValue < 4; jValue++) {
            elementField_host.Flux[jValue][iElement] = 0;//单元数值通量Flux清零，为后面加减做准备
            // f->elements[ie].deltaeig = 0;//每一轮deltaeig清零 针对Roe格式 目前不需要
            // 梯度已经计算，因此不需要再次计算
        }
    }

    // 每条边计算无粘通量，然后根据方向分别加减给两侧单元的Flux。所有边遍历后，所有单元的Flux也就计算出来了
    for (int iedge = 0; iedge < edge_host.num_edge; iedge++) {



        FVM_2D* f = FVM_2D::getInstance();
        int elementL = edge_host.elementL[iedge];
        int elementR = edge_host.elementR[iedge];
        int setID = edge_host.setID[iedge];// 内部edge不属于任何set，因此setID为-1 setID的初始化见readContinueFile
        int boundaryNum = f->boundaryManager.boundaries.size();
        int bType = -1;
        if (setID != -1) {
            // 第iedge个edge的边界类型 = boundaries[setID - 1].type，其中setID = edge_host.setID[iedge]
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
            //对称。对欧拉方程，相当于无粘固壁
        case _BC_symmetry:
            Space::Flux::getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iedge, flux);
            break;
            //无粘固壁
        case _BC_wall_nonViscous:
            Space::Flux::getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iedge, flux);
            break;
            //无滑移绝热
        case _BC_wall_adiabat:
            // TODO: 实现无滑移绝热边界条件的计算
            Space::Flux::getEdgeFlux_wall_adiabat(element_host, elementField_host, edge_host, iedge, flux);
            break;
            //入口
        case _BC_inlet:
            Space::Flux::getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
            break;
            //出口
        case _BC_outlet:
            Space::Flux::getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
            break;
            //远场
        case _BC_inf:
            Space::Flux::getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::inf::ruvp, flux);
            break;
        default:// 内部：bType=-1，边界：bType取_BC_periodic_0到_BC_periodic_9，即6100-6109
            if (elementR != -1) {
                // 周期和内部 统一处理
                Space::Flux::getEdgeFlux_inner_and_periodic(element_host, elementField_host, edge_host, iedge, flux);
            }
        }

        // 更新单元数值通量
        for (int j = 0; j < 4; j++) {
            elementField_host.Flux[j][elementL] += flux[j];
            // 内部边界
            if (bType == -1) {
                // 此处不是elementR==-1，防止周期边界算2次
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
        auto eabs = abs(u * edgeN[0] + v * edgeN[1]);//|euv·en|
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
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getEdgeFlux_wall_adiabat(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, long idx, double* flux) {
    // 固壁边界。函数中idx相当于iedge

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
    // 法向坐标系下 UnL = uL * nx + vL * ny，VnL = -uL * ny + vL * nx
    // UnR = - UnL，VnR = -VnL
    // 原坐标系下 uR = unR*nx-vnR*ny， vR=unR*ny+vnR*nx
    // 因此 uR = -uL 还不如直接取反

    // 按UNITS中BC_WALL.f90的算法，速度取反，可是结果仍然不变
    // 不对 结果云图中没有体现粘性
    //REAL uR = uL - 2 * (uL * nx + vL * ny) * nx;
    //REAL vR = vL - 2 * (uL * nx + vL * ny) * ny;
    REAL uR = -uL;
    REAL vR = -vL;
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    // 计算flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getEdgeFlux_farField(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, long idx, REAL* ruvp_inf, REAL* flux) {
    
    constexpr bool useNew = true;
    if (useNew) {
        getEdgeFlux_farField_2(element, elementField, edge, idx, ruvp_inf, flux);
        return;
    }

    
    // 以下是旧版代码
    REAL nx = edge.normal[0][idx];// normal已初始化
    REAL ny = edge.normal[1][idx];
    REAL x_edge = edge.xy[0][idx];
    REAL y_edge = edge.xy[1][idx];
    REAL length = edge.length[idx];
    long iElementL = edge.elementL[idx];

    // 计算U_L，即edge U的左极限
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementL, U_L);
    // 计算ruvp_L
    REAL ruvp_L[4]{};
    U2NITS::Math::U2ruvp_host(U_L, ruvp_L, GlobalPara::constant::gamma);
    // 检查是否发散
    if (Restrict::outOfRange(ruvp_L)) {
        std::stringstream ss;
        ss << "ElementL=" << iElementL << ", x_edge=" << x_edge << ", y_edge=" << y_edge << ", ";
        ss << "edge=" << idx << ", ";
        ss << "ruvp_L=" << StringProcessor::doubleArray_2_string(ruvp_L, 4) << std::endl;
        ss << "(U2NITS::Space::Flux::getEdgeFlux_farField)\n";
        LogWriter::logAndPrintError(ss.str());
    }

    //根据远场边界条件修正ruvp_L
    Space::Flux::modify_ruvpL_farField(nx, ny, ruvp_L, ruvp_inf);
    U2NITS::Math::ruvp2U_host(ruvp_L, U_L, GlobalPara::constant::gamma);
    // 计算flux
    RiemannSolve(U_L, U_L, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getEdgeFlux_farField_2(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long iEdge, REAL* ruvp_inf, REAL* flux) {
    /*
    2024-04-12 原来的远场边界有问题，无论是激波管还是翼型，远场边界都表现类似固壁边界，即速度为0，仿佛被堵住

    由于激波管在初始化时，内外状态均相同，因此计算得到的Uboundary也跟内状态相同，因此flux=0

    */

    if (!edge_host.has(iEdge)) {
        LogWriter::logAndPrintError("iEdge" + std::to_string(iEdge) + " out of range.\n");
        exit(-1);
    }
    // 远场参数
    real nx = edge_host.normal[0][iEdge];
    real ny = edge_host.normal[1][iEdge];
    real gamma = GlobalPara::constant::gamma;
    real rho_inf = ruvp_inf[0];
    real u_inf = ruvp_inf[1];
    real v_inf = ruvp_inf[2];
    real p_inf = ruvp_inf[3];
    real un_inf = u_inf * nx + v_inf * ny;
    real a2 = gamma * p_inf / rho_inf;
    real Ma2 = (u_inf * u_inf + v_inf * v_inf) / a2;
    real Man2 = (un_inf * un_inf) / a2;
    // 内单元参数
    real UL[4];
    integer element_inner = edge_host.elementL[iEdge];
    for (int i = 0; i < 4; i++) {
        UL[i] = elementField_host.U[i][element_inner];
    }
    // 求边界U
    real Uboundary[4]{};
    bool isGlobalSuperSonic = (Ma2 > 1);
    bool isNormalSuperSonic = (Man2 > 1);
    bool isInlet = (un_inf < 0);
    if (isNormalSuperSonic) {
        if (isInlet) {
            Math::ruvp2U_host(ruvp_inf, Uboundary, gamma);
        }
        else {
            for (int i = 0; i < 4; i++) {
                Uboundary[i] = UL[i];
            }
        }
    }
    else {
        real ruvpL[4]{};
        Math::U2ruvp_host(UL, ruvpL, gamma);
        real rhoa = rho_inf;
        real ua = u_inf;
        real va = v_inf;
        real pa = p_inf;
        real rhod = ruvpL[0];
        real ud = ruvpL[1];
        real vd = ruvpL[2];
        real pd = ruvpL[3];
        real rho0 = rhod;// 参考 CFDPA P264 reference state
        real c02 = gamma * pd / rhod;
        real c0 = sqrt(c02);
        real rho0c0 = rho0 * c0;

        if (isInlet) {
            real pb = 0.5 * (pa + pd - rho0c0 *
                (nx*(ua-ud)+ny*(va-vd)) );
            real pa_pb = pa - pb;
            real rhob = rhoa + (-pa_pb) / c02;
            real ub = ua - nx * pa_pb / rho0c0;
            real vb = va - ny * pa_pb / rho0c0;
            real ruvp_b[4]{ rhob,ub,vb,pb };
            Math::ruvp2U_host(ruvp_b, Uboundary, gamma);
        }
        else {
            real pb = pa;
            real pd_pb = pd - pb;
            real rhob = rhod + (-pd_pb) / c02;
            real ub = ud + nx * pd_pb / rho0c0;
            real vb = vd + ny * (pd_pb) / rho0c0;
            real ruvp_b[4]{ rhob,ub,vb,pb };
            Math::ruvp2U_host(ruvp_b, Uboundary, gamma);
        }
    }
    // 计算flux
    real length = edge_host.length[iEdge];
    RiemannSolve(UL, Uboundary, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, const REAL* ruvp_inf) {
    // 该函数根据远场边界条件修正ruvp_L

    //边界元的Element_R为null，因此nxny一定朝外
    const REAL rho = ruvp[0];
    const REAL u = ruvp[1];
    const REAL v = ruvp[2];
    const REAL p = ruvp[3];
    const REAL rho_inf = ruvp_inf[0];
    const REAL u_inf = ruvp_inf[1];
    const REAL v_inf = ruvp_inf[2];
    const REAL p_inf = ruvp_inf[3];
    const REAL gamma = GlobalPara::constant::gamma;

    //坐标变换
    const REAL rho_n = rho;
    const REAL u_n = u * nx + v * ny;//u cost + v sint
    const REAL v_n = -u * ny + v * nx;//此处有修改。法线向外为正 - u sint + v cost
    const REAL p_n = p;
    const REAL rho_n_inf = rho_inf;
    const REAL p_n_inf = p_inf;
    const REAL u_n_inf = u_inf * nx + v_inf * ny;
    const REAL v_n_inf = -u_inf * ny + v_inf * nx;//此处有修改

    //if (rho_n < 0)system("pause");
    REAL a_n2 = gamma * p_n / rho_n;
    if (a_n2 < 0) {
        std::cout << "Error: a_n2 < 0\n";  
        //system("pause"); 
        //exit(-1); 
    }
    const REAL a_n = sqrt(a_n2);
    REAL a_n_inf2 = gamma * p_n_inf / rho_n_inf;
    if (a_n_inf2 < 0) {
        std::cout << "Error: a_n_inf2 < 0\n";  
        //system("pause");
        //exit(-1); 
    }
    const REAL a_n_inf = sqrt(a_n_inf2);

    REAL v_n_tmp;
    REAL u_n_tmp;
    REAL rho_n_tmp;
    REAL p_n_tmp;

    if (u_n_inf * u_n_inf < a_n_inf * a_n_inf) {//|Vn_inf|<|a_inf|，亚音速
        REAL tmp_a = (gamma - 1.0) / 4.0 * (u_n - u_n_inf) + 0.5 * (a_n + a_n_inf);// 1/2 *( (ga-1)/2 * (u-uinf) + (a+ainf) )
        REAL tmp_a2 = tmp_a * tmp_a;
        if (u_n_inf <= 0.0) {//亚音速入口，提3个边界条件
            REAL tmp_s_inf = p_n_inf / pow(rho_n_inf, gamma);//根据边界rho和p求出熵
            rho_n_tmp = pow(tmp_a2 / tmp_s_inf / gamma, 1.0 / (gamma - 1.0));//边界
            u_n_tmp = 0.5 * (u_n + u_n_inf + 2.0 * a_n / (gamma - 1.0) - 2.0 * a_n_inf / (gamma - 1.0));
            v_n_tmp = v_n_inf;//边界
            p_n_tmp = tmp_a2 * rho_n_tmp / gamma;//边界
        }
        else {//亚音速出口，提1个边界条件 u_n_tmp
            REAL tmp_s = p_n / pow(rho_n, gamma);//内 //pow(x,y):若x<0且y非整数，或者x=0且y<=0，将出现结果错误
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

void U2NITS::Space::Flux::getEdgeFlux_inner_and_periodic(ElementSoA& element, FieldSoA& elementField, EdgeSoA& edge, long idx, REAL* flux) {


    //计算边中点坐标、坐标变换矩阵、坐标逆变换矩阵、边法线方向
    REAL nx = edge.normal[0][idx];// normal已初始化
    REAL ny = edge.normal[1][idx];
    REAL x_edge = edge.xy[0][idx];
    REAL y_edge = edge.xy[1][idx];
    REAL length = edge.length[idx];
    long iElementL = edge.elementL[idx];
    long iElementR = edge.elementR[idx];

    // 计算edge坐标处的U值(分别用两侧单元的分布函数计算)，分为周期边界和内部边界两种情况
    REAL U_L[4]{};
    Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementL, U_L);
    REAL U_R[4]{};
    if (edge.setID[idx] != -1) {
        // 周期边界 应使用对称的edge的坐标计算U值
        int ID = edge.ID[idx];

        //// 法一 用map
        //auto cIter = edge_periodic_pair.find(ID);
        //if (cIter == edge_periodic_pair.end()) {
        //    LogWriter::logAndPrintError("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic\n");
        //    exit(-1);
        //}
        //int ID_pair = cIter->second;

        // 法二 用edge存储
        int ID_pair = edge.periodicPair[ID];
        if (ID_pair < 0 || ID_pair >= edge.num_edge) {
            LogWriter::logAndPrintError("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic\n");
            exit(-1);
        }

        REAL x_edge_pair = edge.xy[0][ID_pair];
        REAL y_edge_pair = edge.xy[1][ID_pair];
        Flux::getUByXYandElementID(element, elementField, x_edge_pair, y_edge_pair, iElementR, U_R);
    }
    else {
        // 内部边界
        Flux::getUByXYandElementID(element, elementField, x_edge, y_edge, iElementR, U_R);
    }


    // 计算flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void U2NITS::Space::Flux::getUByXYandElementID(ElementSoA& element, FieldSoA& elementField, REAL x, REAL y, int elementID, REAL* U_dist) {
    // 该函数依据单元的重构函数，计算某坐标处的U值
    REAL x_elementL = element.xy[0][elementID];
    REAL y_elementL = element.xy[1][elementID];
    // 常量重构
    if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_constant) {
        for (int jValue = 0; jValue < 4; jValue++) {
            U_dist[jValue] = elementField.U[jValue][elementID];
        }
    }
    // 线性重构 根据梯度Ux、Uy计算点(xpoint,ypoint)处_U值
    else if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_linear) {
        for (int jValue = 0; jValue < 4; jValue++) {
            REAL& U_elementL = elementField.U[jValue][elementID];
            REAL& Ux_elementL = elementField.Ux[jValue][elementID];
            REAL& Uy_elementL = elementField.Uy[jValue][elementID];
            U_dist[jValue] = U_elementL + Ux_elementL * (x - x_elementL) + Uy_elementL * (y - y_elementL);
        }
        // 若数据异常，则常量重构
        real ruvp[4]{};
        Math::U2ruvp_host(U_dist, ruvp, GlobalPara::constant::gamma);
        if (Restrict::outOfRange(ruvp)) {
            for (int jValue = 0; jValue < 4; jValue++) {
                U_dist[jValue] = elementField.U[jValue][elementID];
            }
        }
    }
    else {
        LogWriter::logAndPrintError("implemented. @getUByXYandElementID.\n");
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



