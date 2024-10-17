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
该部分代码用于计算通量。
需要用变量edgeField存储各边的数值通量，该变量无需从外部文件读取，而且无需在host与device之间交换。
换言之，对于host计算，只需申请edgeField_host，对于device计算，只需申请edgeField_device
对于host代码，为了避免修改沿途所有函数，最好用全局变量提供访问外部函数的接口
对于device代码，如何直接访问全局函数？
*/
inline bool isUpStreamOfShock_atBoundary_static_inline(myfloat x, myfloat y, myfloat _shock_x, myfloat _shock_y, myfloat _shock_speed, myfloat _t_RK, const myfloat _sqrt3) {
    // CBoundaryDoubleShockReflect的成员函数的静态版本
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
    // 根据单元分布函数，得到某点U值
    for (int i = 0; i < 4; i++) {
        get_U_reconstruct_singleValue(x, y, U_dist[i], iElement, element.xy[0], element.xy[1], elementField.U[i], elementField.Ux[i], elementField.Uy[i], flag_reconstruct);
    }

    // 若数据异常，则常量重构
    myfloat ruvp[4]{};
    U2NITS::Math::U2ruvp_host(U_dist, ruvp, gamma);
    if (U2NITS::Space::Restrict::outOfRange(ruvp)) {
        for (int i = 0; i < 4; i++) {
            U_dist[i] = elementField.U[i][iElement];
        }
    }
}

void get_UR_wallNonViscous(const myfloat* U_L, myfloat* U_R, myfloat nx, myfloat ny) {
    // 滑移壁面，法向速度取反，切向速度相等
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

void get_UR_wall_adiabat(const myfloat* U_L, myfloat* U_R, myfloat nx, myfloat ny) {
    // 无滑移绝热壁面，法向切向速度都取反
    U_R[0] = U_L[0];
    U_R[1] = - U_R[1];
    U_R[2] = - U_R[2];
    U_R[3] = U_L[3];
}

void get_UR_farField(const myfloat* U_L, const myfloat* ruvp_inf, myfloat* U_R, myfloat nx, myfloat ny, myfloat gamma) {
    /*
    参照UNITs和李启兵课件
    以亚声速入口为例，先由内点外推边界点的估计值，
    再建立边界点数值解与估计值的特征相容关系
    解出边界声速，进而得到压力数值解
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
        // -as < qn < 0 注意C++中不能连续写小于号，要拆成 -as<qn && qn<0
        // 亚音速入口
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
    // 检查rho、p是否异常
    if (U2NITS::Space::Restrict::outOfRange(ruvp0)) {
        LogWriter::logAndPrintError("Out of range @get_UR_farField \n");
        CExit::saveAndExit(-1);
    }
    // 根据原始变量ruvp0计算守恒量U_R
    U2NITS::Math::ruvp2U_host(ruvp0, U_R, gamma);
}

void get_UR_inner_and_periodic(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, GPU::EdgeSoA& edge, myint iEdge, myint iElementR, myfloat x_edge, myfloat y_edge, myfloat* U_R, int inviscid_flux_method_flag_reconstruct, myfloat gamma) {
    if (edge.setID[iEdge] != -1) {
        // 周期边界 应使用对称的edge的坐标计算U值
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
        // 内部边界
        get_Uvector_reconstruct_2(element, elementField, x_edge, y_edge, iElementR, U_R, inviscid_flux_method_flag_reconstruct, gamma);
    }

}

void get_flux_RiemannSolve(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny, const myfloat length, myfloat* flux, const int flux_conservation_scheme, myfloat gamma, myfloat rcpcv) {
    // 黎曼求解器，根据UL、UR求解flux
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

inline void clear_element_flux() {
    GPU::GPUSolver2* solver = SolverDataGetter::getSolverInstance();
    GPU::ElementFieldSoA& elementField_host = solver->elementField_host;
    // 单元数值通量Flux清零
    fill_matrix(elementField_host.Flux, 4, elementField_host.num, 0);
}

void edge_convection_flux_kernel(GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, GPU::EdgeFieldSoA& edgeField_host, GPU::BoundarySetMap& boundary_host,
    int flux_conservation_scheme, int inviscid_flux_method_flag_reconstruct, myfloat* inf_ruvp, myfloat* inlet_ruvp, myfloat* outlet_ruvp, myfloat gamma, myfloat rcpcv,
    myfloat dsr_shock_x, myfloat dsr_shock_y, myfloat dsr_shock_speed, myfloat dsr_t_RK, const myfloat _sqrt3, myint iEdge) {
    // 改造成GPU时参数不要传引用
    int setID = edge_host.setID[iEdge];// 内部edge不属于任何set，因此setID为-1 setID的初始化见readContinueFile
    int bType = -1;
    if (setID != -1) {
        // 第iedge个edge的边界类型
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
        // 对称。对欧拉方程，相当于无粘固壁
        get_UR_wallNonViscous(U_L, U_R, nx, ny);
        break;
    case _BC_wall_nonViscous:
        // 无粘固壁
        get_UR_wallNonViscous(U_L, U_R, nx, ny);
        break;
    case _BC_wall_adiabat:
        // 无滑移绝热
        get_UR_wall_adiabat(U_L, U_R, nx, ny);
        break;
    case _BC_inlet:
        // 入口
        get_UR_farField(U_L, inlet_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_outlet:
        // 出口
        get_UR_farField(U_L, outlet_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_inf:
        // 远场
        get_UR_farField(U_L, inf_ruvp, U_R, nx, ny, gamma);
        break;
    case _BC_doubleShockReflect:
        // 双马赫反射
        if (isUpStream_doubleShockReflect) {
            get_UR_farField(U_L, inlet_ruvp, U_R, nx, ny, gamma);
        }
        else {
            get_UR_farField(U_L, outlet_ruvp, U_R, nx, ny, gamma);
        }
        break;

    default:// 内部：bType=-1，边界：bType取_BC_periodic_0到_BC_periodic_9，即6100-6109
        if (iElementR != -1) {
            // 周期和内部 统一处理
            get_UR_inner_and_periodic(element_host, elementField_host, edge_host, iEdge, iElementR, x_edge, y_edge, U_R, inviscid_flux_method_flag_reconstruct, gamma);
        }
    }

    get_flux_RiemannSolve(U_L, U_R, nx, ny, length, flux, flux_conservation_scheme, gamma, rcpcv);

    for (int j = 0; j < 4; j++) {
        edgeField_host.Flux[j][iEdge] = flux[j];
    }

}

// 计算无粘通量 created: 2024-05-01, tgl
void edge_convection_flux() {
    /*
    2023-05-05
    为了改造成GPU代码，需要处理当前代码
    将全局变量变成参数形式
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
    // 以下为双马赫反射的成员函数
    CBoundaryDoubleShockReflect* cbdsr = CBoundaryDoubleShockReflect::getInstance();
    myfloat dsr_shock_x = cbdsr->get_shock_x();
    myfloat dsr_shock_y = cbdsr->get_shock_y();
    myfloat dsr_shock_speed = cbdsr->get_shock_speed();
    myfloat dsr_t_RK = cbdsr->get_t_plus_dt();
    const myfloat _sqrt3 = sqrt(3.0);

    // 每条边计算无粘通量
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

void U2NITS::Space::Flux::calculateFluxHost(GPU::ElementSoA& element_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {
    /*
    20240421: 在保留旧函数calculate_flux_1的同时测试新函数calculate_flux_2，inline避免了函数调用的开销
    函数功能：计算单元数值通量。首先遍历各边，计算边左右守恒量，然后求解黎曼问题得到边通量，最后分别加减到左右单元
    跟之前函数的区别：
    计算得到的Flux不直接加减到邻居单元，而是先存入edgeField_host
    优点1：将大循环拆成一个edge循环和一个element循环，有利于并行
    优点2：将Flux存入edgeField_host，便于插入粘性通量求解部分，为日后加粘性做准备

    先遍历边，计算无粘通量，存进edgeField_host
    然后遍历边，计算粘性通量，对edgeField_host的flux作修改
    然后遍历单元，取edgeField_host的值，加减到自身
    */
    clear_element_flux();// element数值通量清零
    edge_convection_flux();// edge无粘通量
    if (GlobalPara::physicsModel::equation == _EQ_NS) U2NITS::Space::edge_viscous_flux();// edge粘性通量
    add_edge_flux_to_element_flux();// edge数值通量加减到element数值通量

}



