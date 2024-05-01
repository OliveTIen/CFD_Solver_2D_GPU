#include "Flux.h"
#include "../FVM_2D.h"
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
该部分代码用于计算通量。
需要用变量edgeField存储各边的数值通量，该变量无需从外部文件读取，而且无需在host与device之间交换。
换言之，对于host计算，只需申请edgeField_host，对于device计算，只需申请edgeField_device
而且为了避免修改沿途所有函数，最好用全局变量提供访问外部函数的接口
*/

void getUByXYandElementID(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, myfloat x, myfloat y, int elementID, myfloat* U_dist) {
    // 该函数依据单元的重构函数，计算某坐标处的U值
    myfloat x_elementL = element.xy[0][elementID];
    myfloat y_elementL = element.xy[1][elementID];
    // 常量重构
    if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_constant) {
        for (int jValue = 0; jValue < 4; jValue++) {
            U_dist[jValue] = elementField.U[jValue][elementID];
        }
    }
    // 线性重构 根据梯度Ux、Uy计算点(xpoint,ypoint)处_U值
    else if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_linear) {
        for (int jValue = 0; jValue < 4; jValue++) {
            myfloat& U_elementL = elementField.U[jValue][elementID];
            myfloat& Ux_elementL = elementField.Ux[jValue][elementID];
            myfloat& Uy_elementL = elementField.Uy[jValue][elementID];
            U_dist[jValue] = U_elementL + Ux_elementL * (x - x_elementL) + Uy_elementL * (y - y_elementL);
        }
        // 若数据异常，则常量重构
        myfloat ruvp[4]{};

        U2NITS::Math::U2ruvp_host(U_dist, ruvp, GlobalPara::constant::gamma);
        if (U2NITS::Space::Restrict::outOfRange(ruvp)) {
            for (int jValue = 0; jValue < 4; jValue++) {
                U_dist[jValue] = elementField.U[jValue][elementID];
            }
        }
    }
    else {
        LogWriter::logAndPrintError("implemented. @getUByXYandElementID.\n");
    }
}

inline void get_U_reconstruct(myfloat x, myfloat y, myfloat& U_dist, myint i_e, const myfloat* x_e, const myfloat* y_e, const myfloat* U_e, const myfloat* Ux_e, const myfloat* Uy_e, int flag_reconstruct) {
    if (flag_reconstruct == _REC_constant) {
        U_dist = U_e[i_e];
    }
    else if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_linear) {
        U_dist = U_e[i_e] + Ux_e[i_e] * (x - x_e[i_e]) + Uy_e[i_e] * (y - y_e[i_e]);
    }
    else {
        LogWriter::logAndPrintError("unimplemented. @get_U_reconstruct.\n");
    }
}

void get_Uvector_reconstruct_2(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, myfloat x, myfloat y, myint iElement, myfloat* U_dist, int flag_reconstruct) {
    for (int i = 0; i < 4; i++) {
        get_U_reconstruct(x, y, U_dist[i], iElement, element.xy[0], element.xy[1], elementField.U[i], elementField.Ux[i], elementField.Uy[i], flag_reconstruct);
    }

    // 若数据异常，则常量重构
    myfloat ruvp[4]{};
    U2NITS::Math::U2ruvp_host(U_dist, ruvp, GlobalPara::constant::gamma);
    if (U2NITS::Space::Restrict::outOfRange(ruvp)) {
        for (int i = 0; i < 4; i++) {
            U_dist[i] = elementField.U[i][iElement];
        }
    }
}

void RiemannSolve(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny, const myfloat length, myfloat* flux, const int scheme) {
    myfloat faceNormal[2]{ nx,ny };
    myfloat gamma = GlobalPara::constant::gamma;
    switch (scheme) {// GlobalPara::inviscid_flux_method::flux_conservation_scheme
    case _SOL_LocalLaxFriedrichs:
        U2NITS::Space::LocalLaxFriedrichs(UL, UR, nx, ny, length, flux, gamma);
        break;
    case _SOL_Roe:
        U2NITS::Space::ConvectRoeCommon2d(UL, UR, faceNormal, length, flux, gamma);
        break;
    default:
        LogWriter::logAndPrintError("invalid scheme @RiemannSolve");
        exit(-1);
        break;
    }
}

// 滑移壁面，法向速度取反，切向速度相等
void get_UR_wallNonViscous(const myfloat* U_L, myfloat* U_R, myfloat nx, myfloat ny) {
    myfloat uxL = U_L[1] / U_L[0];
    myfloat uyL = U_L[2] / U_L[0];
    myfloat unL = uxL * nx + uyL * ny;// xy为全局坐标系，tn为切向-法向坐标系
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

void getEdgeFlux_wallNonViscous(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, GPU::EdgeSoA& edge, myint iEdge, double* flux) {
    /*
    nx, ny: 归一化法向量
    x_edge, y_edge: 边中点坐标
    length: 边长度
    */

    myfloat nx = edge.normal[0][iEdge];
    myfloat ny = edge.normal[1][iEdge];
    myfloat x_edge = edge.xy[0][iEdge];
    myfloat y_edge = edge.xy[1][iEdge];
    myfloat length = edge.length[iEdge];
    myint iElementL = edge.elementL[iEdge];

    // 计算U_L，即edge U的左极限
    myfloat U_L[4]{};
    get_Uvector_reconstruct_2(element, elementField, x_edge, y_edge, iElementL, U_L, GlobalPara::inviscid_flux_method::flag_reconstruct);
    myfloat U_R[4]{};
    get_UR_wallNonViscous(U_L, U_R, nx, ny);

    // 计算flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

// 无滑移绝热壁面，法向切向速度都取反
void get_UR_wall_adiabat(const myfloat* U_L, myfloat* U_R, myfloat nx, myfloat ny) {
    U_R[0] = U_L[0];
    U_R[1] = - U_R[1];
    U_R[2] = - U_R[2];
    U_R[3] = U_L[3];
}

void getEdgeFlux_wall_adiabat(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, GPU::EdgeSoA& edge, myint iEdge, double* flux) {
    myfloat nx = edge.normal[0][iEdge];
    myfloat ny = edge.normal[1][iEdge];
    myfloat x_edge = edge.xy[0][iEdge];
    myfloat y_edge = edge.xy[1][iEdge];
    myfloat length = edge.length[iEdge];
    myint iElementL = edge.elementL[iEdge];

    myfloat U_L[4]{};
    get_Uvector_reconstruct_2(element, elementField, x_edge, y_edge, iElementL, U_L, GlobalPara::inviscid_flux_method::flag_reconstruct);
    myfloat U_R[4]{};
    get_UR_wall_adiabat(U_L, U_R, nx, ny);

    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void get_UR_farField(const myfloat* U_L, const myfloat* ruvp_inf, myfloat* U_R, myfloat nx, myfloat ny, myfloat gamma) {
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

    // 特征相容关系
    myfloat a2_f = gamma * p_f / rho_f;
    myfloat af = sqrt(a2_f);// 远场声速
    myfloat Ma2_f = (u_f * u_f + v_f * v_f) / a2_f;// 远场马赫数平方
    myfloat ae = sqrt(gamma * p_e / rho_e);// 内点声速
    myfloat qnf = u_f * nx + v_f * ny;// 远场法向速度
    myfloat qne = u_e * nx + v_e * ny;// 内点法向速度
    myfloat rf = qnf - 2.0 * af / ga1;// 第一波特征 R2 = u0n - 2*a0n/(gamma-1)
    myfloat re = qne + 2.0 * ae / ga1;// 第二波特征 u0n - 2*a0n/(gamma-1)
    myfloat qn = 0.5 * (re + rf);
    myfloat as = ga1 * (re - rf) / 4.0;
    myfloat dnt = 0;// 壁面运动速度

    // 判据 q<=-as，超音速入口；-as<q<0，亚音速入口；0<=q<as，亚音速出口；as<=q，超音速出口
    myfloat rho0 = 1, u0 = 0, v0 = 0, p0 = 1;
    if (qn <= -as) {
        // 超音速入口
        rho0 = rho_f;
        u0 = u_f;
        v0 = v_f;
        p0 = p_f;
    }
    else if (qn < 0) {
        // -as < qn < 0 注意C++中不能连续写小于号
        // 亚音速入口 参照UNITs。
        myfloat son = p_f / pow(rho_f, gamma);// 熵
        myfloat qtx = u_f - qnf * nx;
        myfloat qty = v_f - qnf * ny;

        rho0 = pow((as * as / son / gamma), 1.0 / ga1);// 根据S=p/rho^gamma得rho=(p/S)^(1/gamma)=(a2/gamma/S)^(1/gamma)
        u0 = qtx + qn * nx;
        v0 = qty + qn * ny;
        p0 = as * as * rho0 / gamma;// p=rho*R*t=rho/gamma*gamma*R*t=rho/gamma*a2
    }
    else if (qn < as) {
        // 0 <= qn < as
        // 亚音速出口 
        myfloat son = p_e / pow(rho_e, gamma);
        myfloat qtx = u_e - qne * nx;
        myfloat qty = v_e - qne * ny;

        // 后面跟亚音速入口相同
        rho0 = pow((as * as / son / gamma), 1.0 / ga1);
        u0 = qtx + qn * nx;
        v0 = qty + qn * ny;
        p0 = as * as * rho0 / gamma;
    }
    else {
        // as <= q
        // 超音速出口
        rho0 = rho_e;
        u0 = u_e;
        v0 = v_e;
        p0 = p_e;
    }
    myfloat ruvp0[4]{ rho0,u0,v0,p0 };
    // 计算flux
    if (U2NITS::Space::Restrict::outOfRange(ruvp0)) {
        LogWriter::logAndPrintError("Out of range @get_UR_farField \n");
        CExit::saveAndExit(-1);
    }
    U2NITS::Math::ruvp2U_host(ruvp0, U_R, gamma);
}

void getEdgeFlux_farField_3(GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, myint iEdge, myfloat* ruvp_inf, myfloat* flux) {
    /*
    参照UNITs和李启兵课件
    以亚声速入口为例，先由内点外推边界点的估计值，
    再建立边界点数值解与估计值的特征相容关系
    解出边界声速，进而得到压力数值解
    */

    if (!edge_host.has(iEdge)) {
        LogWriter::logAndPrintError("iEdge" + std::to_string(iEdge) + " out of range.\n");
        exit(-1);
    }
    // 远场值
    myfloat nx = edge_host.normal[0][iEdge];
    myfloat ny = edge_host.normal[1][iEdge];
    myfloat gamma = GlobalPara::constant::gamma;
    myfloat ga1 = gamma - 1;
    myfloat rho_f = ruvp_inf[0];
    myfloat u_f = ruvp_inf[1];
    myfloat v_f = ruvp_inf[2];
    myfloat p_f = ruvp_inf[3];

    // 内点值，或边界参数估计值 此处用单元中心值作为估计值
    myfloat U_e[4]{};
    myint element_inner = edge_host.elementL[iEdge];
    for (int i = 0; i < 4; i++) {
        U_e[i] = elementField_host.U[i][element_inner];
    }
    myfloat ruvp_e[4]{};
    U2NITS::Math::U2ruvp_host(U_e, ruvp_e, gamma);
    myfloat rho_e = ruvp_e[0];
    myfloat u_e = ruvp_e[1];
    myfloat v_e = ruvp_e[2];
    myfloat p_e = ruvp_e[3];

    // 特征相容关系
    myfloat a2_f = gamma * p_f / rho_f;
    myfloat af = sqrt(a2_f);// 远场声速
    myfloat Ma2_f = (u_f * u_f + v_f * v_f) / a2_f;// 远场马赫数平方
    myfloat ae = sqrt(gamma * p_e / rho_e);// 内点声速
    myfloat qnf = u_f * nx + v_f * ny;// 远场法向速度
    myfloat qne = u_e * nx + v_e * ny;// 内点法向速度
    myfloat rf = qnf - 2.0 * af / ga1;// 第一波特征 R2 = u0n - 2*a0n/(gamma-1)
    myfloat re = qne + 2.0 * ae / ga1;// 第二波特征 u0n - 2*a0n/(gamma-1)
    myfloat qn = 0.5 * (re + rf);
    myfloat as = ga1 * (re - rf) / 4.0;
    myfloat dnt = 0;// 壁面运动速度

    // 判据 q<=-as，超音速入口；-as<q<0，亚音速入口；0<=q<as，亚音速出口；as<=q，超音速出口
    myfloat rho0 = 1, u0 = 0, v0 = 0, p0 = 1;
    if (qn <= -as) {
        // 超音速入口
        rho0 = rho_f;
        u0 = u_f;
        v0 = v_f;
        p0 = p_f;
    }
    else if (qn < 0) {
        // -as < qn < 0 注意C++中不能连续写小于号
        // 亚音速入口 参照UNITs。
        myfloat son = p_f / pow(rho_f, gamma);// 熵
        myfloat qtx = u_f - qnf * nx;
        myfloat qty = v_f - qnf * ny;

        rho0 = pow((as * as / son / gamma), 1.0 / ga1);// 根据S=p/rho^gamma得rho=(p/S)^(1/gamma)=(a2/gamma/S)^(1/gamma)
        u0 = qtx + qn * nx;
        v0 = qty + qn * ny;
        p0 = as * as * rho0 / gamma;// p=rho*R*t=rho/gamma*gamma*R*t=rho/gamma*a2
    }
    else if (qn < as) {
        // 0 <= qn < as
        // 亚音速出口 
        myfloat son = p_e / pow(rho_e, gamma);
        myfloat qtx = u_e - qne * nx;
        myfloat qty = v_e - qne * ny;

        // 后面跟亚音速入口相同
        rho0 = pow((as * as / son / gamma), 1.0 / ga1);
        u0 = qtx + qn * nx;
        v0 = qty + qn * ny;
        p0 = as * as * rho0 / gamma;
    }
    else {
        // as <= q
        // 超音速出口
        rho0 = rho_e;
        u0 = u_e;
        v0 = v_e;
        p0 = p_e;
    }
    myfloat ruvp0[4]{ rho0,u0,v0,p0 };
    // 计算flux
    myfloat U0[4]{};
    U2NITS::Math::ruvp2U_host(ruvp0, U0, gamma);
    if (U2NITS::Space::Restrict::outOfRange(ruvp0)) {
        LogWriter::logAndPrintError("Out of range @U2NITS::Space::Flux::getEdgeFlux_farField_3 \n");
        CExit::saveAndExit(-1);
    }

    myfloat length = edge_host.length[iEdge];
    RiemannSolve(U_e, U0, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);

}

void getEdgeFlux_farField_4(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, GPU::EdgeSoA& edge, myint iEdge, myfloat* ruvp_inf, myfloat* flux) {
    /*
    参照UNITs和李启兵课件
    以亚声速入口为例，先由内点外推边界点的估计值，
    再建立边界点数值解与估计值的特征相容关系
    解出边界声速，进而得到压力数值解

    20240421: 修改代码，添加get_UR_farField，使其与其他代码形式统一
    */
    myfloat nx = edge.normal[0][iEdge];
    myfloat ny = edge.normal[1][iEdge];
    myfloat x_edge = edge.xy[0][iEdge];
    myfloat y_edge = edge.xy[1][iEdge];
    myfloat length = edge.length[iEdge];
    myint iElementL = edge.elementL[iEdge];

    // 计算U_L，即edge U的左极限
    myfloat U_L[4]{};
    get_Uvector_reconstruct_2(element, elementField, x_edge, y_edge, iElementL, U_L, _REC_constant);
    myfloat U_R[4]{};
    get_UR_farField(U_L, ruvp_inf, U_R, nx, ny, GlobalPara::constant::gamma);

    // 计算flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void getEdgeFlux_farField(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, GPU::EdgeSoA& edge, myint iEdge, myfloat* ruvp_inf, myfloat* flux) {
    getEdgeFlux_farField_4(element, elementField, edge, iEdge, ruvp_inf, flux);
}

void getEdgeFlux_farField_doubleShockReflect(GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, myint iEdge, myfloat* flux) {
    /*
    计算边界元是位于上游还是下游，并分别用入口和出口的参数作为远场边界参数
    */
    myfloat edgex = edge_host.xy[0][iEdge];
    myfloat edgey = edge_host.xy[1][iEdge];
    CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
    bool isUpStream = pCDS->isUpStreamOfShock_atBoundary(edgex, edgey, pCDS->get_t_plus_dt());
    if (isUpStream) {
        getEdgeFlux_farField(element_host, elementField_host, edge_host, iEdge,
            GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
    }
    else {
        getEdgeFlux_farField(element_host, elementField_host, edge_host, iEdge,
            GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
    }
}



void getEdgeFlux_inner_and_periodic(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, GPU::EdgeSoA& edge, long idx, myfloat* flux) {


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
    get_Uvector_reconstruct_2(element, elementField, x_edge, y_edge, iElementL, U_L, GlobalPara::inviscid_flux_method::flag_reconstruct);
    myfloat U_R[4]{};
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
            LogWriter::logAndPrintError("periodicPairNotFoundException, @GPU::Space::Flux::getEdgeFlux_inner_and_periodic_kernel\n");
            exit(-1);
        }

        myfloat x_edge_pair = edge.xy[0][ID_pair];
        myfloat y_edge_pair = edge.xy[1][ID_pair];
        get_Uvector_reconstruct_2(element, elementField, x_edge_pair, y_edge_pair, iElementR, U_R, GlobalPara::inviscid_flux_method::flag_reconstruct);
    }
    else {
        // 内部边界
        get_Uvector_reconstruct_2(element, elementField, x_edge, y_edge, iElementR, U_R, GlobalPara::inviscid_flux_method::flag_reconstruct);
    }


    // 计算flux
    RiemannSolve(U_L, U_R, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

inline void fill_matrix(myfloat* mat[], int nrow, myint ncol, myfloat _Val) {
    /*
    输入：mat通常为4xn矩阵。数据类型为myfloat* mat[4]，4个指针各自指向一个一维数组，因此各行非连续分配，要分别memset
    n较大(可能超过21亿)，因此用myint

    memset以字节为单位进行赋值，所赋值的范围是0x00～0xFF: memset(mat[i], 0, ncol * sizeof(myfloat));
    https://zhuanlan.zhihu.com/p/479101175
    int赋值为0：参数_Val给0即可
    int赋值为1：参数_Val不能直接给1，因为实际上表示的整数为0x01010101。没有比较好的方法
    float赋值为0：

    还是推荐用std::fill，能赋任意值，比memset安全，比loop快
    https://iq.opengenus.org/set-array-to-0-in-cpp/
    */
    for (int i = 0; i < nrow; i++) {
        std::fill(mat[i], mat[i] + ncol, _Val);
    }
}

inline void clear_edge_flux() {
    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementFieldSoA& elementField_host = solver->elementField_host;
    // 单元数值通量Flux清零
    fill_matrix(elementField_host.Flux, 4, elementField_host.num, 0);
}

// 旧。先计算通量，然后分别加减到相邻单元
inline void calculate_flux_1(GPU::ElementSoA& element_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {
    /*
    参考：void Solver_2D::calFlux()
    TODO:
    计算周期边界的代码需要修改。最多支持10对周期边界，若有更多，会按inner计算，导致错误
    */

    // 单元数值通量Flux清零，为后面加减做准备. 梯度已经计算，因此不需要再次计算
    fill_matrix(elementField_host.Flux, 4, elementField_host.num, 0);

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

        switch (bType) {
        case _BC_symmetry:
            // 对称。对欧拉方程，相当于无粘固壁
            getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iedge, flux);
            break;
        case _BC_wall_nonViscous:
            // 无粘固壁
            getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iedge, flux);
            break;
        case _BC_wall_adiabat:
            // 无滑移绝热壁
            getEdgeFlux_wall_adiabat(element_host, elementField_host, edge_host, iedge, flux);
            break;
        case _BC_inlet:
            // 入口
            getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
            break;
        case _BC_outlet:
            // 出口
            getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
            break;
        case _BC_inf:
            // 远场
            getEdgeFlux_farField(element_host, elementField_host, edge_host, iedge,
                GlobalPara::boundaryCondition::_2D::inf::ruvp, flux);
            break;
        case _BC_doubleShockReflect:
            // 双马赫反射
            getEdgeFlux_farField_doubleShockReflect(element_host, elementField_host, edge_host, iedge, flux);

            break;
        default:// 内部：bType=-1，边界：bType取_BC_periodic_0到_BC_periodic_9，即6100-6109
            if (elementR != -1) {
                // 周期和内部 统一处理
                getEdgeFlux_inner_and_periodic(element_host, elementField_host, edge_host, iedge, flux);
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

/*
计算无粘通量
created: 2024-05-01, tgl
*/
void edge_convection_flux() {
    if (GlobalPara::physicsModel::equation) {

    }

    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementSoA& element_host = solver->element_host;
    GPU::ElementFieldSoA& elementField_host = solver->elementField_host;
    GPU::EdgeSoA& edge_host = solver->edge_host;
    GPU::EdgeFieldSoA& edgeField_host = solver->edgeField_host;
    // 每条边计算无粘通量
    for (int iEdge = 0; iEdge < edge_host.num_edge; iEdge++) {

        FVM_2D* f = FVM_2D::getInstance();
        int elementL = edge_host.elementL[iEdge];
        int elementR = edge_host.elementR[iEdge];
        int setID = edge_host.setID[iEdge];// 内部edge不属于任何set，因此setID为-1 setID的初始化见readContinueFile
        int boundaryNum = f->boundaryManager.boundaries.size();
        int bType = -1;
        if (setID != -1) {
            // 第iedge个edge的边界类型 = boundaries[setID - 1].type，其中setID = edge_host.setID[iedge]
            bType = f->boundaryManager.boundaries[setID - 1].type;
        }

        myfloat nx = edge_host.normal[0][iEdge];
        myfloat ny = edge_host.normal[1][iEdge];
        myfloat x_edge = edge_host.xy[0][iEdge];
        myfloat y_edge = edge_host.xy[1][iEdge];
        myfloat length = edge_host.length[iEdge];
        myint iElementL = edge_host.elementL[iEdge];

        myfloat flux[4]{};

        switch (bType) {
        case _BC_symmetry:
            // 对称。对欧拉方程，相当于无粘固壁
            getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iEdge, flux);
            break;
        case _BC_wall_nonViscous:
            // 无粘固壁
            getEdgeFlux_wallNonViscous(element_host, elementField_host, edge_host, iEdge, flux);
            break;
        case _BC_wall_adiabat:
            // 无滑移绝热
            getEdgeFlux_wall_adiabat(element_host, elementField_host, edge_host, iEdge, flux);
            break;
        case _BC_inlet:
            // 入口
            getEdgeFlux_farField(element_host, elementField_host, edge_host, iEdge,
                GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
            break;
        case _BC_outlet:
            // 出口
            getEdgeFlux_farField(element_host, elementField_host, edge_host, iEdge,
                GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
            break;
        case _BC_inf:
            // 远场
            getEdgeFlux_farField(element_host, elementField_host, edge_host, iEdge,
                GlobalPara::boundaryCondition::_2D::inf::ruvp, flux);
            break;
        case _BC_doubleShockReflect:
            // 双马赫反射
            getEdgeFlux_farField_doubleShockReflect(element_host, elementField_host, edge_host, iEdge, flux);

            break;
        default:// 内部：bType=-1，边界：bType取_BC_periodic_0到_BC_periodic_9，即6100-6109
            if (elementR != -1) {
                // 周期和内部 统一处理
                getEdgeFlux_inner_and_periodic(element_host, elementField_host, edge_host, iEdge, flux);
            }
        }

        for (int j = 0; j < 4; j++) {
            edgeField_host.Flux[j][iEdge] = flux[j];
        }
    }// end for iEdge

}

inline void add_edge_flux_to_element_flux() {
    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementSoA& element_host = solver->element_host;
    GPU::ElementFieldSoA& elementField_host = solver->elementField_host;
    GPU::EdgeSoA& edge_host = solver->edge_host;
    GPU::EdgeFieldSoA& edgeField_host = solver->edgeField_host;

    // 加减到相邻单元
    for (myint iElement = 0; iElement < element_host.num_element; iElement++) {

        // 参见Gradient
        // 求前单元体积、edgeID、edge朝向正负、外单元ID
        myfloat volumeC = element_host.volume[iElement];// 体积
        // 初始化edgeID 超出数组范围的非法值赋为-1
        int edgeID[4]{ -1,-1,-1,-1 };// 第iElement个element的第i条边的ID
        for (int i = 0; i < 4; i++) {
            edgeID[i] = element_host.edges[i][iElement];
            if (!edge_host.has(edgeID[i])) {
                edgeID[i] = -1;
            }
        }
        // 初始化边朝向edgeSign
        int edgeSign[4]{ 0,0,0,0 };// 1表示边朝外，-1表示边朝内
        for (int i = 0; i < 4; i++) {
            int edgeIDi = edgeID[i];
            if (edgeIDi == -1) {
                continue;
            }
            if (iElement != edge_host.elementR[edgeIDi]) {
                edgeSign[i] = 1;// currentElement=elementL，朝外
            }
            else {
                edgeSign[i] = -1;// currentElement=elementR，朝内
            }
        }
        // 根据边朝向，加减到中心单元。向外为加，向内为减
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

/*
新函数，计算通量。
先遍历边，计算无粘通量，存进edgeField_host
然后遍历边，计算粘性通量，对edgeField_host的flux作修改
然后遍历单元，取edgeField_host的值，加减到自身
*/
inline void calculate_flux_2() {
    /*
    20240421: 在保留旧函数calculate_flux_1的同时测试新函数calculate_flux_2，inline避免了函数调用的开销
    函数功能：计算单元数值通量。首先遍历各边，计算边左右守恒量，然后求解黎曼问题得到边通量，最后分别加减到左右单元
    跟之前函数的区别：
    计算得到的Flux不直接加减到邻居单元，而是先存入edgeField_host
    优点1：将大循环拆成一个edge循环和一个element循环，有利于并行
    优点2：将Flux存入edgeField_host，便于插入粘性通量求解部分，为日后加粘性做准备
    */
    clear_edge_flux();// edge数值通量清零
    edge_convection_flux();// edge无粘通量
    if (GlobalPara::physicsModel::equation == _EQ_NS) U2NITS::Space::edge_viscous_flux();// edge粘性通量
    add_edge_flux_to_element_flux();// edge数值通量加减到element数值通量

}

void U2NITS::Space::Flux::calculateFluxHost(GPU::ElementSoA& element_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {
    calculate_flux_2();
}



