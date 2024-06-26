#include "Solver_2D.h"
#include "FVM_2D.h"
#include "output/LogWriter.h"
#include "space/RiemannSolver.h"
#include "space/Reconstructor.h"

//myfloat Solver_2D::RK3alpha[6]{ 0.0,1.0 / 4.0,1.0 / 6.0,3.0 / 8.0 ,0.5,1.0 };
//myfloat Solver_2D::RK5alpha[6]{ 0.0,1.0 / 4.0,1.0 / 6.0,3.0 / 8.0 ,0.5,1.0 };


//Eigen::Matrix4d Solver_2D::get_matrixT(myfloat nx, myfloat ny, int flag) {
//    Eigen::Matrix4d matT;
//    ny *= flag;//逆矩阵就是sin变为相反数
//    matT <<
//        1, 0, 0, 0,
//        0, nx, ny, 0,
//        0, -ny, nx, 0,
//        0, 0, 0, 1;
//    return matT;
//}

//void Solver_2D::U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, myfloat& lambda) {
//    myfloat rho = U[0];
//    myfloat u = U[1] / U[0];
//    myfloat v = U[2] / U[0];
//    myfloat E = U[3] / U[0];
//    myfloat gamma = GlobalPara::constant::gamma;
//    Math_2D::U_2_F(U, F, gamma);
//    myfloat p = Math_2D::get_p(rho, gamma, E, u, v);
//    lambda = sqrt(u * u + v * v) + sqrt(gamma * p / rho);
//
//}

void Solver_2D::evolve(myfloat dt) {
    if (GlobalPara::time::time_advance == _EVO_explicit) {
        evolve_explicit(dt);
    }
    else if (GlobalPara::time::time_advance == _EVO_rk3) {
        evolve_RK3(dt);
    }
    else {
        LogWriter::logAndPrintError("Invalid input \"time.time_advance\". Will exit.");
    }
}

void Solver_2D::evolve_explicit(myfloat dt) {
    //数值通量
    calFlux();
    //时间推进
    FVM_2D* f = FVM_2D::getInstance();
    myfloat omega;//面积
//#pragma omp parallel for
    for (int ie = 0; ie < f->elements.size(); ie++) {
        omega = f->elements[ie].calArea(f);
        for (int j = 0; j < 4; j++) {
            f->elements[ie].U[j] -= dt / omega * f->elements[ie].Flux[j];
        }
    }
}

void Solver_2D::evolve_RK3(myfloat dt) {

    FVM_2D* f = FVM_2D::getInstance();
    //定义(局部)结构体，用于存储U值
    struct Data {
        myfloat U[4]{};
    };
    std::vector<Data> elements_old(f->elements.size());
    std::vector<Data> k1(elements_old), k2(k1), k3(k1);
    //复制
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            elements_old[ie].U[j] = f->elements[ie].U[j];
        }
    }
    
    //计算面积并存储，防止反复计算
    std::vector<myfloat>areas; areas.resize(f->elements.size());
    for (int ie = 0; ie < f->elements.size(); ie++) {
        areas[ie] = f->elements[ie].calArea(f);
    }

    
    //1.k1=f(t_n,y_n)
    //1.1.基于U值计算数值通量
    calFlux();
    //1.2.求k1
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            k1[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / areas[ie];//注意负号
        }
    }
    //2.k2=f(t_n+dt/2,y_n+dt/2*k1)
    //2.1.求U值y_n+dt/2*k1
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].U[j] = elements_old[ie].U[j] + dt/2.0 * k1[ie].U[j];
        }
    }
    //2.2.基于U值计算数值通量
    calFlux();
    //2.3.求k2
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            k2[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / areas[ie];//注意负号
        }
    }
    //3.k3=f(t_n+dt,y_n-dt*k1+2*dt*k2)
    //3.1.求U值y_n-dt*k1+2*dt*k2
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].U[j] = elements_old[ie].U[j] - dt * k1[ie].U[j] + 2.0 * dt * k2[ie].U[j];
        }
    }
    //3.2.基于U值计算数值通量
    calFlux();
    //3.3.求k3
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            k3[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / areas[ie];//注意负号
        }
    }
    //4.y_n+1=y_n+dt/6.0*(k1+4*k2+k3)
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].U[j] = elements_old[ie].U[j] + dt / 6.0 * (k1[ie].U[j] + 4 * k2[ie].U[j] + k3[ie].U[j]);
        }
    }


}


void Solver_2D::calFlux() {
    /*
    已初始化的参数：
        单元坐标

    计算通量过程：
	1. 初始化所有单元数值通量为0
    2. 计算单元梯度
        根据单元坐标、单元U、邻居坐标、邻居U，计算梯度
        限制器修正梯度
        检查异常值
        输入：单元坐标、单元U、邻居坐标、邻居U
        输出：梯度
    3. 计算边界数值通量
        计算边法向量、边坐标
        计算边左右单元在边界的U值(左右极限UL、UR)
        根据边类型，修正UL、UR
        计算黎曼问题
            输入：UL UR nx ny edgeLength
            输出：边界的数值通量
    4. 计算单元数值通量
        将边界数值通量加到单元数值通量上
    5. 显式时间推进

    */


    FVM_2D* f = FVM_2D::getInstance();

//#pragma omp parallel for
    //准备工作。初始化单元数值通量、分布函数
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].Flux[j] = 0;//单元数值通量Flux清零，为后面加减做准备
            f->elements[ie].deltaeig = 0;//每一轮deltaeig清零
            if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_constant
                && GlobalPara::physicsModel::equation == _EQ_euler) {
                // 如果是常量重构+Euler，则无需计算梯度
            }
            else {
                // 否则需要计算梯度，且用限制器修正
                Reconstructor::Element_T3_updateSlope_Barth(f, &(f->elements[ie]));
                //f->elements[ie].updateSlope_Barth(f);//线性重构，Barth限制器
            }
        }
    }

//#pragma omp parallel for
    //每条边计算无粘通量，然后根据方向分别加减给两侧单元的Flux。所有边遍历后，所有单元的Flux也就计算出来了
    for (int iedge = 0; iedge < f->edges.size(); iedge++) {
        ////计算每条边的黎曼通量
        myfloat flux[4]{};

        {
            FVM_2D* f = FVM_2D::getInstance();
            Edge_2D* pE = &(f->edges[iedge]);
            if (pE->pElement_R == nullptr) {//边界
                const int bType = f->boundaryManager.boundaries[pE->setID - 1].type;
                switch (bType) {
                    //对称。对欧拉方程而言，对称边界即无粘固壁边界
                case _BC_symmetry:
                    getEdgeFlux_wallNonViscous(pE, flux);
                    break;

                    //无粘固壁
                case _BC_wall_nonViscous:
                    getEdgeFlux_wallNonViscous(pE, flux);
                    break;
                    //无滑移绝热
                case _BC_wall_adiabat:
                    
                    break;
                    //入口
                case _BC_inlet:
                    getEdgeFlux_farfield(pE, GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
                    break;
                    //出口
                case _BC_outlet:
                    getEdgeFlux_farfield(pE, GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
                    break;
                    //远场
                case _BC_inf:
                    getEdgeFlux_farfield(pE, GlobalPara::boundaryCondition::_2D::inf::ruvp, flux);
                    break;
                }
                //周期
                if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) getEdgeFlux_periodic(pE, flux);//最多10对周期边界
            }
            else//内部
            {
                getEdgeFlux_inner(pE, flux);
            }

        }

        ////减给两侧单元的Flux
        for (int j = 0; j < 4; j++) {
            if (f->edges[iedge].pElement_R != nullptr) {//内部edge
                f->edges[iedge].pElement_L->Flux[j] += flux[j];
                f->edges[iedge].pElement_R->Flux[j] -= flux[j];
            }
            else {//边界edge
                f->edges[iedge].pElement_L->Flux[j] += flux[j];
            }
        }
    }

}

void Solver_2D::getEdgeFlux_inner(Edge_2D* pE, myfloat* flux) {
    //功能：计算内部边界的数值通量
    
    myfloat x_edge, y_edge;     
    myfloat U_L[4], U_R[4];
    myfloat nx, ny;
    pE->getxy(FVM_2D::getInstance(), x_edge, y_edge);
    pE->pElement_L->get_U(x_edge, y_edge, U_L);
    pE->pElement_R->get_U(x_edge, y_edge, U_R);
    pE->getDirectionN(nx, ny);
    RiemannSolve(U_L, U_R, nx, ny, pE->getLength(), flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);

}

void Solver_2D::getEdgeFlux_wallNonViscous(Edge_2D* pE, myfloat* flux) {

    myfloat x_edge, y_edge;
    myfloat U_L[4];
    myfloat nx, ny;
    pE->getDirectionN(nx, ny);
    pE->getxy(FVM_2D::getInstance(), x_edge, y_edge);
    pE->pElement_L->get_U(x_edge, y_edge, U_L);

    //计算U_R
    myfloat uL = U_L[1] / U_L[0];
    myfloat vL = U_L[2] / U_L[0];
    myfloat uR = uL - 2 * (uL * nx + vL * ny) * nx;
    myfloat vR = vL - 2 * (uL * nx + vL * ny) * ny;
    myfloat U_R[4];
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    RiemannSolve(U_L, U_R, nx, ny, pE->getLength(), flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

void Solver_2D::getEdgeFlux_farfield(Edge_2D* pE, const myfloat* ruvp_inf, myfloat* flux) {
    //函数功能：根据远场边界条件计算边界数值通量
    FVM_2D* f = FVM_2D::getInstance();
    myfloat nx, ny; 
    myfloat x_edge, y_edge;
    myfloat U_L[4];
    myfloat ruvp_L[4];
    //计算ruvp_L
    pE->getDirectionN(nx, ny);
    pE->getxy(FVM_2D::getInstance(), x_edge, y_edge);
    pE->pElement_L->get_U(x_edge, y_edge, U_L);
    Math_2D::U_2_ruvp(U_L, ruvp_L, GlobalPara::constant::gamma);
    //检查ruvp_L合法性 检测是否发散
    if (ruvp_L[3] < 0) {
        std::cout << "ElementID=" << pE->pElement_L->ID << ",x=" << pE->pElement_L->x << ",y=" << pE->pElement_L->y
            << ":\tp<0" << "\t(Edge_2D::calFlux_Riemann)\n";
    }
    if (ruvp_L[0] < 0) {
        std::cout << "ElementID=" << pE->pElement_L->ID << ",x=" << pE->pElement_L->x << ",y=" << pE->pElement_L->y
            << ":\trho<0" << "\t(Edge_2D::calFlux_Riemann)\n";
    }
    if (ruvp_L[0] > 900) {
        std::cout << "Warning rho>900, ElementID=" << pE->pElement_L->ID << "\n";
    }
    //根据远场边界条件修正ruvp_L
    modify_ruvpL_in_farfield(nx, ny, ruvp_L, ruvp_inf);
    Math_2D::ruvp_2_U(ruvp_L, U_L, GlobalPara::constant::gamma);
    RiemannSolve(U_L, U_L, nx, ny, pE->getLength(), flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);

}

void Solver_2D::getEdgeFlux_periodic(Edge_2D* pE, myfloat* flux) {
    //FVM_2D* f = FVM_2D::getInstance();
    ////计算边中点坐标、坐标变换矩阵、坐标逆变换矩阵、边法线方向
    //myfloat x_edge, y_edge; pE->getxy(f, x_edge, y_edge);
    //myfloat nx, ny;    pE->getDirectionN(nx, ny);
    ////计算U_L、U_R
    //Eigen::Vector4d U_L = pE->pElement_L->get_U(x_edge, y_edge);
    //Eigen::Vector4d U_R;
    //{
    //    //找到这条edge对应的虚拟pElement_R(本来该edge没有pElement_R，但我们根据周期边界找到它的虚拟pElement_R)
    //    Edge_2D* pE_pair = f->boundaryManager.get_pairEdge_periodic(pE);//找到edge在周期边界中对应的edge_1
    //    Element_2D* pElement_R = pE_pair->pElement_L;//pEdge_1的pElement_L即为pEdge的虚拟pElement_R
    //    //获取edge_1的中点坐标
    //    myfloat x_virtual_edge, y_virtual_edge;
    //    pE_pair->getxy(f, x_virtual_edge, y_virtual_edge);
    //    //重构。根据虚拟单元计算U_R
    //    U_R = pElement_R->get_U(x_virtual_edge, y_virtual_edge);
    //}
    ////计算数值通量向量。其元素为各守恒量的数值通量
    //myfloat UL[4]{ U_L[0],U_L[1],U_L[2],U_L[3] };
    //myfloat UR[4]{ U_R[0],U_R[1],U_R[2],U_R[3] };
    //RiemannSolve(UL, UR, nx, ny, pE->getLength(), flux,
    //    GlobalPara::inviscid_flux_method::flux_conservation_scheme);

    // 重新写
    FVM_2D* f = FVM_2D::getInstance();
    myfloat x_edge{};
    myfloat y_edge{};
    pE->getxy(f, x_edge, y_edge);
    myfloat nx{};
    myfloat ny{};
    pE->getDirectionN(nx, ny);
    myfloat UL[4]{};
    pE->pElement_L->get_U(x_edge, y_edge, UL);
    Edge_2D* pE_pair = f->boundaryManager.get_pairEdge_periodic(pE);
    Element_2D* pElement_pairL = pE_pair->pElement_L;
    myfloat x_edge_pair{};
    myfloat y_edge_pair{};
    pE_pair->getxy(f, x_edge_pair, y_edge_pair);
    myfloat UR[4]{};
    pElement_pairL->get_U(x_edge_pair,y_edge_pair,UR);
    RiemannSolve(UL, UR, nx, ny, pE->getLength(), flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
}

//void Solver_2D::calFlux_Roe_2() {
//    //功能：计算单元数值通量
//
//    //准备工作
//    FVM_2D* f = FVM_2D::getInstance();
//    for (int ie = 0; ie < f->elements.size(); ie++) {
//        for (int j = 0; j < 4; j++) {
//            f->elements[ie].Flux[j] = 0;//单元数值通量Flux清零，为后面加减做准备
//            if (GlobalStatic::flag_reconstruct == _REC_linear)f->elements[ie].updateSlope_Barth(f);//线性重构，即更新单元分布函数
//        }
//    }
//
//    //单元deltaeig清零，然后赋值
//    Compute_Deltaeig();
//
//    for (int ie = 0; ie < f->edges.size(); ie++) {
//        Edge_2D* pE = &(f->edges[ie]);
//
//        //2023-8-16 新方法
//        //myfloat NL = Neighbor(1, i);左右单元
//        //myfloat NR = Neighbor(2, i);
//        //myfloat DX = VECX(1, i);
//        //myfloat DY = VECX(2, i);
//        //myfloat n1 = FNode(1, i);节点坐标
//        //myfloat n2 = FNode(2, i);
//
//
//        //myfloat xmid, ymid;//边中点坐标
//        //pE->getxy(f, xmid, ymid);
//        //myfloat sav1n, sav2n;//边法向？？？？？？？？？？？问问老师
//        ////pE->getDirectionN(sav1n, sav2n);
//        //myfloat sav1 = (f->getNodeByID(pE->nodes[1])->x - f->getNodeByID(pE->nodes[0])->x);
//        //myfloat sav2 = (f->getNodeByID(pE->nodes[1])->y - f->getNodeByID(pE->nodes[0])->y);
//        //myfloat sav = sqrt(sav1 * sav1 + sav2 * sav2);
//        //sav1n = sav1 / sav;
//        //sav2n = sav2 / sav;
//
//        ////新增 直接导致发散
//        myfloat sav1, sav2, sav1n, sav2n, sav;
//        pE->getDirectionN(sav1, sav2);
//        sav = sqrt(sav1 * sav1 + sav2 * sav2);
//        sav1n = sav1 / sav;
//        sav2n = sav2 / sav;
//
//
//        //do ig=1, NG 高斯积分点
//
//        //计算左ruvp
//        //!计算左单元在边中点处的ruvp
//        Eigen::Vector4d UL = pE->get_UL();
//        myfloat _UL[4]{ UL[0],UL[1],UL[2],UL[3] };
//        myfloat ruvp_L[4];
//        Math_2D::U_2_ruvp(_UL, ruvp_L, GlobalPara::constant::gamma);
//        myfloat rl = ruvp_L[0];
//        myfloat ul = ruvp_L[1];
//        myfloat vl = ruvp_L[2];
//        myfloat pl = ruvp_L[3];
//        //计算左单元中心点处的ruvp
//        myfloat* UL_ele = pE->pElement_L->U;
//        myfloat ruvp_L_ele[4];
//        Math_2D::U_2_ruvp(UL_ele, ruvp_L_ele, GlobalPara::constant::gamma);
//        myfloat rlt = ruvp_L_ele[0];
//        myfloat ult = ruvp_L_ele[1];
//        myfloat vlt = ruvp_L_ele[2];
//        myfloat plt = ruvp_L_ele[3];
//        //限制器。若边界处的值太离谱，则用单元中心处的值替代
//        if (rl < 0 || pl < 0 || rl>100000 || pl>1000000) {
//            rl = rlt;
//            ul = ult;
//            vl = vlt;
//            pl = plt;
//        }
//        myfloat d_small = 1.0e-16;
//        rl = max(rl, d_small);
//        pl = max(pl, d_small);
//        rlt = max(rlt, d_small);
//        plt = max(plt, d_small);
//
//        //计算右ruvp
//        myfloat rr;
//        myfloat ur;
//        myfloat vr;
//        myfloat pr;
//        myfloat rrt;
//        myfloat urt;
//        myfloat vrt;
//        myfloat prt;
//        //判断边界类型
//        if (pE->pElement_R != nullptr) {//内部
//
//            //!计算右单元在边中点处的ruvp
//            Eigen::Vector4d UR = pE->get_UR();
//            myfloat _UR[4]{ UR[0],UR[1],UR[2],UR[3] };
//            myfloat ruvp_R[4];
//            Math_2D::U_2_ruvp(_UR, ruvp_R, GlobalPara::constant::gamma);
//            rr = ruvp_R[0];
//            ur = ruvp_R[1];
//            vr = ruvp_R[2];
//            pr = ruvp_R[3];
//            //计算右单元中心点处的ruvp
//            myfloat* UR_ele = pE->pElement_R->U;
//            myfloat ruvp_R_ele[4];
//            Math_2D::U_2_ruvp(UR_ele, ruvp_R_ele, GlobalPara::constant::gamma);
//            rrt = ruvp_R_ele[0];
//            urt = ruvp_R_ele[1];
//            vrt = ruvp_R_ele[2];
//            prt = ruvp_R_ele[3];
//            //限制器
//            if (rr < 0 || pr < 0 || rr>100000 || pr>1000000) {
//                rr = rrt;
//                ur = urt;
//                vr = vrt;
//                pr = prt;
//            }
//        }
//
//        else {//边界 尚未完成
//            const int bType = f->boundaryManager.boundaries[pE->setID - 1].type;
//            switch (bType) {
//                //对称。对欧拉方程而言，对称边界即无粘固壁边界
//            case _BC_symmetry:
//                //无需break
//
//                //无粘固壁
//            case _BC_wall_nonViscous:
//            {
//
//            }
//            //无滑移绝热
//            case _BC_wall_adiabat:
//                
//                //入口
//            case _BC_inlet:
//                //return calEdgeFlux_LLF_farfield(pE, GlobalPara::boundaryCondition::_2D::inlet::ruvp);
//                //出口
//            case _BC_outlet:
//                //return calEdgeFlux_LLF_farfield(pE, GlobalPara::boundaryCondition::_2D::outlet::ruvp);
//                //远场
//            case _BC_inf:
//                //return calEdgeFlux_LLF_farfield(pE, GlobalPara::boundaryCondition::_2D::inf::ruvp);
//                std::cout << "Error: 未完成此部分代码(Solver_2D::calFlux_Roe_2)\n";
//                break;
//
//            default:
//                break;
//            }
//            //周期
//            if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {//最多10对周期边界
//                //计算边中点坐标、坐标变换矩阵、坐标逆变换矩阵、边法线方向
//                myfloat x_edge, y_edge; pE->getxy(f, x_edge, y_edge);
//                std::vector<myfloat> directionN = pE->getDirectionN(f);
//                myfloat nx = directionN[0];
//                myfloat ny = directionN[1];
//                Eigen::Matrix4d T = get_matrixT(nx, ny, 1);
//                Eigen::Matrix4d T_inv = get_matrixT(nx, ny, -1);
//                //计算U_L、U_R
//                Eigen::Vector4d U_L = pE->pElement_L->get_U(x_edge, y_edge);
//                Eigen::Vector4d U_R;
//                myfloat* UR_ele;
//                {
//                    //找到这条edge对应的虚拟pElement_R(本来该edge没有pElement_R，但我们根据周期边界找到它的虚拟pElement_R)
//                    Edge_2D* pEdge_1 = f->boundaryManager.get_pairEdge_periodic(pE);//找到edge在周期边界中对应的edge_1
//                    Element_2D* virtual_pElement_R = pEdge_1->pElement_L;//pEdge_1的pElement_L即为pEdge的虚拟pElement_R
//                    UR_ele = virtual_pElement_R->U;
//                    //获取edge_1的中点坐标
//                    myfloat x_virtual_edge, y_virtual_edge;
//                    pEdge_1->getxy(f, x_virtual_edge, y_virtual_edge);
//                    //重构。根据虚拟单元计算U_R
//                    U_R = virtual_pElement_R->get_U(x_virtual_edge, y_virtual_edge);
//                }
//
//                //计算右单元在边中点处的ruvp
//                Eigen::Vector4d UR = U_R;
//                myfloat _UR[4]{ UR[0],UR[1],UR[2],UR[3] };
//                myfloat ruvp_R[4];
//                Math_2D::U_2_ruvp(_UR, ruvp_R, GlobalPara::constant::gamma);
//                rr = ruvp_R[0];
//                ur = ruvp_R[1];
//                vr = ruvp_R[2];
//                pr = ruvp_R[3];
//                //计算右单元中心点处的ruvp
//                myfloat ruvp_R_ele[4];
//                Math_2D::U_2_ruvp(UR_ele, ruvp_R_ele, GlobalPara::constant::gamma);
//                rrt = ruvp_R_ele[0];
//                urt = ruvp_R_ele[1];
//                vrt = ruvp_R_ele[2];
//                prt = ruvp_R_ele[3];
//                //限制器
//                if (rr < 0 || pr < 0 || rr>100000 || pr>1000000) {
//                    rr = rrt;
//                    ur = urt;
//                    vr = vrt;
//                    pr = prt;
//                }
//
//
//            }
//        }
//
//        rr = max(rr, d_small);
//        pr = max(pr, d_small);
//        rrt = max(rrt, d_small);
//        prt = max(prt, d_small);
//
//        const myfloat& gama = GlobalPara::constant::gamma;
//        const myfloat ga1 = gama - 1;
//        myfloat vaml = ult * ult + vlt * vlt;
//        myfloat vamr = urt * urt + vrt * vrt;
//        myfloat hlt = plt * gama / (rlt * (gama - 1.)) + 0.5 * vaml;
//        myfloat hrt = prt * gama / (rrt * (gama - 1.)) + 0.5 * vamr;
//        myfloat acl = sqrt(gama * plt / rlt);
//        myfloat acr = sqrt(gama * prt / rrt);
//
//        ////Roe平均
//        myfloat rrorl = sqrt(rrt / rlt);
//        myfloat rrorlp1 = 1. + rrorl;
//        myfloat rm = sqrt(rrt * rlt);
//        myfloat um = (ult + urt * rrorl) / rrorlp1;
//        myfloat vm = (vlt + vrt * rrorl) / rrorlp1;
//        myfloat hm = (hlt + hrt * rrorl) / rrorlp1;
//        myfloat vm2 = um * um + vm * vm;
//        myfloat am2 = ga1 * abs(hm - 0.5 * vm2);
//        myfloat am = sqrt(am2);
//        myfloat am2i = 1. / am2;
//        myfloat ami = 1. / am;
//        //特征值
//
//        myfloat unormal = sav1 * um + sav2 * vm;
//        myfloat anormal = am * sav;
//        myfloat dacou = acr - acl;
//        myfloat eig1 = abs(unormal);
//        myfloat eig2 = abs(unormal + anormal);  //, deigen)
//        myfloat eig3 = abs(unormal - anormal); //, deigen)
//
//
//        ////entropy fix, H
//        myfloat yita = pE->pElement_L->deltaeig;//使用Compute_Deltaeig()进行初始化，该函数未完成
//        if (pE->pElement_R != nullptr) yita = max(yita, pE->pElement_R->deltaeig);
//        myfloat epss = 2. * yita * sav;
//
//        //!!Harten entropy fix
//        if (eig1 < epss) eig1 = 0.5 * (eig1 * eig1 / epss + epss);
//        if (eig2 < epss) eig2 = 0.5 * (eig2 * eig2 / epss + epss);
//        if (eig3 < epss) eig3 = 0.5 * (eig3 * eig3 / epss + epss);
//
//        //计算hl hr
//        vaml = ul * ul + vl * vl;
//        vamr = ur * ur + vr * vr;
//        myfloat hl = pl * gama / (rl * (gama - 1.)) + 0.5 * vaml;
//        myfloat hr = pr * gama / (rr * (gama - 1.)) + 0.5 * vamr;
//
//        //inviscid flux
//        myfloat dunormal = (sav1n * (ur - ul) + sav2n * (vr - vl));
//        myfloat alph1 = eig1 * (rr - rl - (pr - pl) * am2i);
//        myfloat alph2 = 0.5 * am2i * eig2 * (pr - pl + rm * am * dunormal);
//        myfloat alph3 = 0.5 * am2i * eig3 * (pr - pl - rm * am * dunormal);
//        myfloat alph4 = alph1 + alph2 + alph3;
//        myfloat alph5 = am * (alph2 - alph3);
//        myfloat alph6 = eig1 * rm * (ur - ul - sav1n * dunormal);
//        myfloat alph7 = eig1 * rm * (vr - vl - sav2n * dunormal);
//
//        myfloat ulnormal = sav1 * ul + sav2 * vl;
//        myfloat urnormal = sav1 * ur + sav2 * vr;
//        myfloat rulnormal = rl * ulnormal;
//        myfloat rurnormal = rr * urnormal;
//        myfloat unormaln = (unormal) / sav;
//
//        //数值通量
//        myfloat flux[4];
//        flux[0] = 0.5 * (rulnormal + rurnormal - alph4);
//        flux[1] = 0.5 * (rulnormal * ul + rurnormal * ur
//            + sav1 * (pl + pr) - um * alph4 - sav1n * alph5 - alph6);
//        flux[2] = 0.5 * (rulnormal * vl + rurnormal * vr
//            + sav2 * (pl + pr) - vm * alph4 - sav2n * alph5 - alph7);
//        flux[3] = 0.5 * (rulnormal * hl + rurnormal * hr
//            - hm * alph4 - unormaln * alph5
//            - um * alph6 - vm * alph7 + am2 * alph1 / ga1);
//        
//        //将边的数值通量加减到单元上
//        for (int j = 0; j < 4; j++) {
//            pE->pElement_L->Flux[j] += flux[j];
//            if (pE->pElement_R != nullptr) {//内部edge
//                pE->pElement_R->Flux[j] -= flux[j];
//            }
//        }
//
//
//    }
//
//    //debug
//    for (int ie = 0; ie < f->edges.size(); ie++) {
//        if (f->edges[ie].ID == 5374) {
//            std::string str = GlobalStatic::myfloatArray_2_string(f->edges[ie].pElement_L->U, 4) + "\n";
//            GlobalStatic::writeLogAndCout(str);
//        }
//    }
//
//}

void Solver_2D::Compute_Deltaeig() {
    //准备工作
    FVM_2D* f = FVM_2D::getInstance();
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].deltaeig = 0;//每一轮deltaeig清零
        }
    }
    //计算每个单元的deltaeig
    const myfloat gamma = GlobalPara::constant::gamma;
    for (int ie = 0; ie < f->edges.size(); ie++) {
        Edge_2D* pE = &(f->edges[ie]);

        //myfloat xmid, ymid;//边中点坐标
        //pE->getxy(f, xmid, ymid);
        //myfloat sav1n, sav2n;//边法向？？？？？？？？？？？问问老师
        ////pE->getDirectionN(sav1n, sav2n);
        //myfloat sav1 = (f->getNodeByID(pE->nodes[1])->x - f->getNodeByID(pE->nodes[0])->x);
        //myfloat sav2 = (f->getNodeByID(pE->nodes[1])->y - f->getNodeByID(pE->nodes[0])->y);
        //myfloat sav = sqrt(sav1 * sav1 + sav2 * sav2);
        //sav1n = sav1 / sav;
        //sav2n = sav2 / sav;

        ////新增 直接导致发散
        myfloat sav1, sav2, sav1n, sav2n, sav;
        pE->getDirectionN(sav1, sav2);
        sav = sqrt(sav1 * sav1 + sav2 * sav2);
        sav1n = sav1 / sav;
        sav2n = sav2 / sav;


        //计算左ruvp
        //计算左单元在边中点处的ruvp
        myfloat rl, ul, vl, pl, acl;
        myfloat rr, ur, vr, pr, acr;
        if (pE->pElement_R != nullptr) {
            //计算单元在边中点处的ruvp
            myfloat _UL[4]{};
            myfloat _UR[4]{};
            pE->get_ULUR(_UL, _UR);
            myfloat ruvp_L[4];
            myfloat ruvp_R[4];
            Math_2D::U_2_ruvp(_UL, ruvp_L, GlobalPara::constant::gamma);
            Math_2D::U_2_ruvp(_UR, ruvp_R, GlobalPara::constant::gamma);
            rl = ruvp_L[0];
            ul = ruvp_L[1];
            vl = ruvp_L[2];
            pl = ruvp_L[3];
            acl = sqrt(gamma * pl / rl);
            rr = ruvp_R[0];
            ur = ruvp_R[1];
            vr = ruvp_R[2];
            pr = ruvp_R[3];
            acr = sqrt(gamma * pr / rr);
            //计算deltaeig
            myfloat deig = abs(sav1n * (ul - ur) + sav2n * (vl - vr)) + abs(acl - acr);
            pE->pElement_L->deltaeig = max(pE->pElement_L->deltaeig, deig);
            pE->pElement_R->deltaeig = max(pE->pElement_R->deltaeig, deig);
        }
    }
}

void Solver_2D::RiemannSolve(const myfloat* UL, const myfloat* UR, const myfloat nx, 
    const myfloat ny, const myfloat length, myfloat* flux, const int conservation_scheme) {
    int returnCode = RiemannSolver::solve(UL, UR, nx, ny, length, flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);
    if (returnCode == RiemannSolver::ReturnStatus::invalid_solver_type) {
        LogWriter::logAndPrintError("invalid RiemannSolver type.\n Will exit.\n");
        exit(returnCode);
    }
    else if (returnCode == RiemannSolver::ReturnStatus::compute_error) {
        LogWriter::logAndPrintError("compute error in RiemannSolver.\n Will exit.\n");
        exit(returnCode);
    }
}

void Solver_2D::modify_ruvpL_in_farfield(const myfloat nx, const myfloat ny, myfloat* ruvp, const myfloat* ruvp_inf) {
    //边界元的Element_R为null，因此nxny一定朝外
    const myfloat rho = ruvp[0];
    const myfloat u = ruvp[1];
    const myfloat v = ruvp[2];
    const myfloat p = ruvp[3];
    const myfloat rho_inf = ruvp_inf[0];
    const myfloat u_inf = ruvp_inf[1];
    const myfloat v_inf = ruvp_inf[2];
    const myfloat p_inf = ruvp_inf[3];
    const myfloat gamma = GlobalPara::constant::gamma;

    //坐标变换
    const myfloat rho_n = rho;
    const myfloat u_n = u * nx + v * ny;//u cost + v sint
    const myfloat v_n = -u * ny + v * nx;//此处有修改。法线向外为正 - u sint + v cost
    //const myfloat v_n = u * ny - v * nx;
    const myfloat p_n = p;
    const myfloat rho_n_inf = rho_inf;
    const myfloat p_n_inf = p_inf;
    const myfloat u_n_inf = u_inf * nx + v_inf * ny;
    const myfloat v_n_inf = -u_inf * ny + v_inf * nx;//此处有修改

    //if (rho_n < 0)system("pause");
    myfloat a_n2 = gamma * p_n / rho_n;
    if (a_n2 < 0) { std::cout << "Error: a_n2 < 0\n";  system("pause"); }
    const myfloat a_n = sqrt(a_n2);
    myfloat a_n_inf2 = gamma * p_n_inf / rho_n_inf;
    if (a_n_inf2 < 0) { std::cout << "Error: a_n_inf2 < 0\n";  system("pause"); }
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

void Solver_2D::flux_viscous(Edge_2D* pE, myfloat* flux) {
    // 应放进Solver_2D::calFlux()的for each face boundary中
    // 此处，不同边界条件的处理全放进一个函数中
    // edge又称为face



    FVM_2D* pFVM2D = FVM_2D::getInstance();
    myfloat nx, ny;
    pE->getDirectionN(nx, ny);
    myfloat x_edge, y_edge;
    pE->getxy(pFVM2D, x_edge, y_edge);
    const myfloat& gamma = GlobalPara::constant::gamma;// 引用，避免占空间
    const int nVar = 4;
    myfloat sav = pE->getRefLength();// ! 此处已经修改！前面的梯度张量部分也要修改！
    myfloat sav1 = sav * nx;
    myfloat sav2 = sav * ny;
    myfloat& sav1n = nx;
    myfloat& sav2n = ny;
    myfloat& dx = sav1;
    myfloat& dy = sav2;
    bool isWall = false;

    //// UL
    myfloat U_L[4]{};
    myfloat ruvp_L[4]{};
    myfloat rl{};// 初始化为0.0
    myfloat ul{};// myfloat rl(ruvp_L[0])相当于 myfloat rl = ruvp_L[0]
    myfloat vl{};
    myfloat pl{};
    pE->pElement_L->get_U(x_edge, y_edge, U_L);
    Math_2D::U_2_ruvp(U_L, ruvp_L, gamma);
    rl = ruvp_L[0];
    ul = ruvp_L[1];
    vl = ruvp_L[2];
    pl = ruvp_L[3];

    //// UR
    myfloat U_R[4]{};
    myfloat ruvp_R[4]{};
    myfloat rr{};// 初始化为0.0
    myfloat ur{};
    myfloat vr{};
    myfloat pr{};
    // inner face
    if (pE->pElement_R != nullptr) {
        pE->pElement_R->get_U(x_edge, y_edge, U_R);
        Math_2D::U_2_ruvp(U_R, ruvp_R, gamma);
        rr = ruvp_R[0];
        ur = ruvp_R[1];
        vr = ruvp_R[2];
        pr = ruvp_R[3];

    }
    // boundary face
    else {
        const int bType = pFVM2D->boundaryManager.boundaries[pE->setID - 1].type;

        // 固壁
        if (bType == _BC_wall_adiabat ||
            bType == _BC_wall_isothermal ||
            bType == _BC_wall_nonViscous) {
            ul = 0;
            vl = 0;
            rr = rl;
            ur = ul;
            vr = vl;
            pr = pl;
            isWall = true;
        }
        else {
            rr = rl;
            ur = ul;
            vr = vl;
            pr = pl;
        }
    }

    //// viscous flux
    myfloat r = 0.5 * (rl + rr);
    myfloat u = 0.5 * (ul + ur);
    myfloat v = 0.5 * (vl + vr);
    myfloat p = 0.5 * (pl + pr);
    // gradient
    myfloat gradC[2][nVar]{};
    myfloat tgd[2][nVar]{};
    myfloat yita = 0.5;
    
    for (int k = 0; k < nVar; k++) {
        gradC[0][k] = pE->pElement_L->Ux[k];
        gradC[1][k] = pE->pElement_L->Uy[k];
        tgd[0][k] = pE->pElement_L->Ux[k] * dx;
        tgd[1][k] = pE->pElement_L->Uy[k] * dy;
    }
    if (pE->pElement_R != nullptr) {
        myfloat dxk = (std::min)(pE->pElement_L->calArea(pFVM2D), pE->pElement_R->calArea(pFVM2D)) / sav;
        for (int k = 0; k < nVar; k++) {
            // !!!! 注意WangQ的"flux_viscous.f90第142行"不是用的ruvp是用的puvt，要订正，因为后面耗散是t的耗散
            tgd[0][k] = 0.5 * (tgd[0][k] + pE->pElement_R->Ux[k] * dx)
                + yita * sav1n * (ruvp_R[k] - ruvp_L[k]) / dxk;
            tgd[1][k] = 0.5 * (tgd[1][k] + pE->pElement_R->Uy[k] * dy)
                + yita * sav2n * (ruvp_R[k] - ruvp_L[k]) / dxk;
        }
    }
    // 接下来是WangQ的"flux_viscous.f90第150行
    myfloat drdx(tgd[0][0]);
    myfloat dudx(tgd[0][1]);
    myfloat dvdx(tgd[0][2]);
    myfloat dpdx(tgd[0][3]);
    myfloat drdy(tgd[1][0]);
    myfloat dudy(tgd[1][1]);
    myfloat dvdy(tgd[1][2]);
    myfloat dpdy(tgd[1][3]);
    if (isWall) {

    }

}

void Solver_2D::flux_viscous_2(Edge_2D* pE, myfloat* flux_viscous) {
    FVM_2D* pFVM2D = FVM_2D::getInstance();
    bool isWall = false;
    float reflen = pE->getRefLength();
    myfloat x_edge{}, y_edge{};
    pE->getxy(pFVM2D, x_edge, y_edge);

    myfloat vvh[2]{ (x_edge - pE->pElement_L->x) / reflen ,
        (y_edge - pE->pElement_L->y) / reflen };
    // reflen 无量纲化的长度 取的是绝对值
    //// UL
    myfloat U_L[4]{}, U_R[4]{};
    myfloat ruvp_L[4]{}, ruvp_R[4]{};
    const myfloat gamma = GlobalPara::constant::gamma;
    const myfloat R = GlobalPara::constant::R;
    pE->pElement_L->get_U(x_edge, y_edge, U_L);
    Math_2D::U_2_ruvp(U_L, ruvp_L, gamma);
    myfloat rl = ruvp_L[0];
    myfloat ul = ruvp_L[1];
    myfloat vl = ruvp_L[2];
    myfloat pl = ruvp_L[3];
    myfloat tl = pl / rl / R; // p=rho*R*T

    if (rl < 0 || pl < 0 || tl < 0) {
        // do nothing
    }
    myfloat rr{}, ur{}, vr{}, pr{}, tr{};
    if (pE->pElement_R != nullptr) {
        pE->pElement_R->get_U(x_edge, y_edge, U_R);
        Math_2D::U_2_ruvp(U_R, ruvp_R, gamma);
        rr = ruvp_R[0];
        ur = ruvp_R[1];
        vr = ruvp_R[2];
        pr = ruvp_R[3];
        tr = pr / rr / R;
    }
    else {
        const int bType = pFVM2D->
            boundaryManager.boundaries[pE->setID - 1].type;

        // 固壁
        if (bType == _BC_wall_adiabat ||
            bType == _BC_wall_isothermal ||
            bType == _BC_wall_nonViscous) {
            rr = rl;
            pr = pl;
            tr = tl;
            ul = 0;
            vl = 0;
            ur = 0;
            vr = 0;
            isWall = true;
        }
        else {
            rr = rl;
            pr = pl;
            tr = tl;

            ur = ul;
            vr = vl;
        }
    }

    // ------------
    // viscous flux
    // ------------
    myfloat tem = 0.5 * (tl + tr);
    myfloat uu = 0.5 * (ul + ur);
    myfloat vv = 0.5 * (vl + vr);
    myfloat pp = 0.5 * (pl + pr);
    rr = pp / (R * tem);// p=rho*R*T
    // Cp-Cv=R, Cp/Cv=gamma
    // p=rho*R*T,R=287.06
    // air: Cp=1005 [J/(kg*K)] 
    //      Cv=717.94
    // 经验证，air的Cp Cv满足以上两个关系式

    // gradient
    const int nVar = 4;
    myfloat gradC[2][nVar]{};
    myfloat tgd[2][nVar]{};
    myfloat yita = 0.5;

    for (int k = 0; k < nVar; k++) {
        gradC[0][k] = pE->pElement_L->Ux[k];
        gradC[1][k] = pE->pElement_L->Uy[k];
        //tgd[0][k] = pE->pElement_L->Ux[k] * dx;
        //tgd[1][k] = pE->pElement_L->Uy[k] * dy;
    }
}
