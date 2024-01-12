#include "Solver_2D.h"
#include "FVM_2D.h"

//double Solver_2D::RK3alpha[6]{ 0.0,1.0 / 4.0,1.0 / 6.0,3.0 / 8.0 ,0.5,1.0 };
//double Solver_2D::RK5alpha[6]{ 0.0,1.0 / 4.0,1.0 / 6.0,3.0 / 8.0 ,0.5,1.0 };


//Eigen::Matrix4d Solver_2D::get_matrixT(double nx, double ny, int flag) {
//    Eigen::Matrix4d matT;
//    ny *= flag;//逆矩阵就是sin变为相反数
//    matT <<
//        1, 0, 0, 0,
//        0, nx, ny, 0,
//        0, -ny, nx, 0,
//        0, 0, 0, 1;
//    return matT;
//}

void Solver_2D::U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, double& lambda) {
    double rho = U[0];
    double u = U[1] / U[0];
    double v = U[2] / U[0];
    double E = U[3] / U[0];
    double gamma = Constant::gamma;
    Math_2D::U_2_F(U, F, gamma);
    double p = Math_2D::get_p(rho, gamma, E, u, v);
    lambda = sqrt(u * u + v * v) + sqrt(gamma * p / rho);

}

void Solver_2D::evolve(double dt) {
    evolve_explicit(dt);
    //evolve_RK3(dt);//平均一个文件38.5s 每步385/531.0s
}

void Solver_2D::evolve_explicit(double dt) {
    //数值通量
    calFlux();
    //时间推进
    FVM_2D* f = FVM_2D::pFVM2D;
    double omega;//面积
    for (int ie = 0; ie < f->elements.size(); ie++) {
        omega = f->elements[ie].calArea(f);
        for (int j = 0; j < 4; j++) {
            f->elements[ie].U[j] -= dt / omega * f->elements[ie].Flux[j];
        }
    }
}

void Solver_2D::evolve_RK3(double dt) {

    FVM_2D* f = FVM_2D::pFVM2D;
    //定义(局部)结构体，用于存储U值
    struct Data {
        double U[4]{};
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
    std::vector<double>omegas; omegas.resize(f->elements.size());
    for (int ie = 0; ie < f->elements.size(); ie++) {
        omegas[ie] = f->elements[ie].calArea(f);
    }

    
    //1.k1=f(t_n,y_n)
    //1.1.基于U值计算数值通量
    calFlux();
    //1.2.求k1
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            k1[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / omegas[ie];//注意负号
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
            k2[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / omegas[ie];//注意负号
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
            k3[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / omegas[ie];//注意负号
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

    calFlux_current();
    return;
}

void Solver_2D::calFlux_current() {
    FVM_2D* f = FVM_2D::pFVM2D;

    //准备工作。初始化单元数值通量、分布函数
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].Flux[j] = 0;//单元数值通量Flux清零，为后面加减做准备
            f->elements[ie].deltaeig = 0;//每一轮deltaeig清零
            if (global::flag_reconstruct == _REC_linear)f->elements[ie].updateSlope_Barth(f);//线性重构，即更新单元分布函数
        }
    }

    //每条边计算黎曼通量，然后根据方向分别加减给两侧单元的Flux。所有边遍历后，所有单元的Flux也就计算出来了
    for (int iedge = 0; iedge < f->edges.size(); iedge++) {
        ////计算每条边的黎曼通量
        double flux[4];

        {
            FVM_2D* f = FVM_2D::pFVM2D;
            Edge_2D* pE = &(f->edges[iedge]);
            if (pE->pElement_R == nullptr) {//边界
                const int bType = f->boundaryManager.vBoundarySets[pE->setID - 1].type;
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

void Solver_2D::getEdgeFlux_inner(Edge_2D* pE, double* flux) {
    //功能：计算内部边界的数值通量
    
    double x_edge, y_edge;     
    double U_L[4], U_R[4];
    double nx, ny;
    pE->getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    pE->pElement_L->get_U(x_edge, y_edge, U_L);
    pE->pElement_R->get_U(x_edge, y_edge, U_R);
    pE->getDirectionN(nx, ny);
    LLF_new_current(U_L, U_R, nx, ny, pE->getLength(), flux);

}

void Solver_2D::getEdgeFlux_wallNonViscous(Edge_2D* pE, double* flux) {

    double x_edge, y_edge;
    double U_L[4];
    double nx, ny;
    pE->getDirectionN(nx, ny);
    pE->getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    pE->pElement_L->get_U(x_edge, y_edge, U_L);

    //计算U_R
    double uL = U_L[1] / U_L[0];
    double vL = U_L[2] / U_L[0];
    double uR = uL - 2 * (uL * nx + vL * ny) * nx;
    double vR = vL - 2 * (uL * nx + vL * ny) * ny;
    double U_R[4];
    U_R[0] = U_L[0];
    U_R[1] = U_R[0] * uR;
    U_R[2] = U_R[0] * vR;
    U_R[3] = U_L[3];

    LLF_new_current(U_L, U_R, nx, ny, pE->getLength(), flux);
}

void Solver_2D::getEdgeFlux_farfield(Edge_2D* pE, const double* ruvp_inf, double* flux) {
    //函数功能：根据远场边界条件计算边界数值通量
    FVM_2D* f = FVM_2D::pFVM2D;
    double nx, ny; 
    double x_edge, y_edge;
    double U_L[4];
    double ruvp_L[4];
    //计算ruvp_L
    pE->getDirectionN(nx, ny);
    pE->getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    pE->pElement_L->get_U(x_edge, y_edge, U_L);
    Math_2D::U_2_ruvp(U_L, ruvp_L, Constant::gamma);
    //检查ruvp_L合法性
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
    cal_ruvp_farfield_new(nx, ny, ruvp_L, ruvp_inf);
    Math_2D::ruvp_2_U(ruvp_L, U_L, Constant::gamma);
    LLF_new_current(U_L, U_L, nx, ny, pE->getLength(), flux);

}

void Solver_2D::getEdgeFlux_periodic(Edge_2D* pE, double* flux) {
    FVM_2D* f = FVM_2D::pFVM2D;
    //计算边中点坐标、坐标变换矩阵、坐标逆变换矩阵、边法线方向
    double x_edge, y_edge; pE->getxy(f, x_edge, y_edge);
    double nx, ny;    pE->getDirectionN(nx, ny);
    //计算U_L、U_R
    Eigen::Vector4d U_L = pE->pElement_L->get_U(x_edge, y_edge);
    Eigen::Vector4d U_R;
    {
        //找到这条edge对应的虚拟pElement_R(本来该edge没有pElement_R，但我们根据周期边界找到它的虚拟pElement_R)
        Edge_2D* pEdge_1 = f->boundaryManager.get_pairEdge_periodic(pE);//找到edge在周期边界中对应的edge_1
        Element_T3* virtual_pElement_R = pEdge_1->pElement_L;//pEdge_1的pElement_L即为pEdge的虚拟pElement_R
        //获取edge_1的中点坐标
        double x_virtual_edge, y_virtual_edge;
        pEdge_1->getxy(f, x_virtual_edge, y_virtual_edge);
        //重构。根据虚拟单元计算U_R
        U_R = virtual_pElement_R->get_U(x_virtual_edge, y_virtual_edge);
    }
    //计算数值通量向量。其元素为各守恒量的数值通量
    double UL[4]{ U_L[0],U_L[1],U_L[2],U_L[3] };
    double UR[4]{ U_R[0],U_R[1],U_R[2],U_R[3] };
    LLF_new_current(UL, UR, nx, ny, pE->getLength(), flux);
}

void Solver_2D::LLF_new_current(const double* UL, const double* UR, const double nx, const double ny, const double length, double* flux) {
    //功能：根据ULUR等参数，计算flux

    const double& gamma = Constant::gamma;
    double ruvpL[4], ruvpR[4];
    Math_2D::U_2_ruvp(UL, ruvpL, gamma);
    Math_2D::U_2_ruvp(UR, ruvpR, gamma);
    double rL = ruvpL[0];
    double uL = ruvpL[1];
    double vL = ruvpL[2];
    double pL = ruvpL[3];
    double rR = ruvpR[0];
    double uR = ruvpR[1];
    double vR = ruvpR[2];
    double pR = ruvpR[3];
    double vaml = uL * uL + vL * vL;
    double vamr = uR * uR + vR * vR;
    double unL = nx * uL + ny * vL;
    double unR = nx * uR + ny * vR;

    double FnL[4], FnR[4];
    FnL[0] = rL * unL;
    FnL[1] = rL * unL * uL + nx * pL;
    FnL[2] = rL * unL * vL + ny * pL;
    double hl = pL * gamma / (rL * (gamma - 1.)) + 0.5 * vaml;
    FnL[3] = rL * unL * hl;
    FnR[0] = rR * unR;
    FnR[1] = rR * unR * uR + nx * pR;
    FnR[2] = rR * unR * vR + ny * pR;
    double hr = pR * gamma / (rR * (gamma - 1.)) + 0.5 * vamr;
    FnR[3] = rR * unR * hr;
    double lambdaMax = max(abs(unL) + sqrt(gamma * pL / rL), abs(unR) + sqrt(gamma * pR / rR));

    for (int i = 0; i < 4; i++) {
        flux[i] = 0.5 * (FnL[i] + FnR[i] - lambdaMax * (UR[i] - UL[i]));
        flux[i] *= length;
    }


    
    //////测试成功 new_1_old 老师的fortran代码
    //const double& gamma = Constant::gamma;
    //double nx_len, ny_len, sav1n, sav2n, len;
    //sav1n = nx;
    //sav2n = ny;
    //len = length;
    //nx_len = length * nx;
    //ny_len = length * ny;

    //double ruvpL[4], ruvpR[4];
    //Math_2D::U_2_ruvp(UL, ruvpL, gamma);
    //Math_2D::U_2_ruvp(UR, ruvpR, gamma);
    //double rL = ruvpL[0];
    //double uL = ruvpL[1];
    //double vL = ruvpL[2];
    //double pL = ruvpL[3];
    //double rR = ruvpR[0];
    //double uR = ruvpR[1];
    //double vR = ruvpR[2];
    //double pR = ruvpR[3];
    //
    //const double ga1 = gamma - 1;
    //double vaml = uL * uL + vL * vL;
    //double vamr = uR * uR + vR * vR;
    //double hl = pL * gamma / (rL * (gamma - 1.)) + 0.5 * vaml;
    //double hr = pR * gamma / (rR * (gamma - 1.)) + 0.5 * vamr;
    //double el = pL / (rL * ga1) + 0.5 * vaml;
    //double er = pR / (rR * ga1) + 0.5 * vamr;
    //double acl = sqrt(gamma * pL / rL);
    //double acr = sqrt(gamma * pR / rR);

    //double unL_len = (nx_len * uL + ny_len * vL);
    //double unR_len = (nx_len * uR + ny_len * vR);

    //double coex = max(abs(unL_len), abs(unR_len)) + max(acl, acr) * len;//lambdaMax=max(lambdaL,lambdaR)

    //double runL_len = rL * unL_len;
    //double runR_len = rR * unR_len;

    //flux[0] = 0.5 * (runL_len + runR_len - coex * (rR - rL));
    //flux[1] = 0.5 * (runL_len * uL + runR_len * uR + nx_len * (pL + pR) - coex * (rR * uR - rL * uL));
    //flux[2] = 0.5 * (runL_len * vL + runR_len * vR + ny_len * (pL + pR) - coex * (rR * vR - rL * vL));
    //flux[3] = 0.5 * (runL_len * hl + runR_len * hr - coex * (rR * er - rL * el));
}

void Solver_2D::LLF_test_2(const double* UL, const double* UR, const double nx, const double ny, const double length, double* flux) {
    //功能：根据ULUR等参数，计算flux

    const double& gamma = Constant::gamma;
    double ruvpL[4], ruvpR[4];
    Math_2D::U_2_ruvp(UL, ruvpL, gamma);
    Math_2D::U_2_ruvp(UR, ruvpR, gamma);
    double rL = ruvpL[0];
    double uL = ruvpL[1];
    double vL = ruvpL[2];
    double pL = ruvpL[3];
    double rR = ruvpR[0];
    double uR = ruvpR[1];
    double vR = ruvpR[2];
    double pR = ruvpR[3];
    double vaml = uL * uL + vL * vL;
    double vamr = uR * uR + vR * vR;
    double unL = nx * uL + ny * vL;
    double unR = nx * uR + ny * vR;

    double FnL[4], FnR[4];
    FnL[0] = rL * unL;
    FnL[1] = rL * unL * uL + nx * pL;
    FnL[2] = rL * unL * vL + ny * pL;
    double hl = pL * gamma / (rL * (gamma - 1.)) + 0.5 * vaml;
    FnL[3] = rL * unL * hl;
    FnR[0] = rR * unR;
    FnR[1] = rR * unR * uR + nx * pR;
    FnR[2] = rR * unR * vR + ny * pR;
    double hr = pR * gamma / (rR * (gamma - 1.)) + 0.5 * vamr;
    FnR[3] = rR * unR * hr;
    double lambdaMax = max(abs(unL) + sqrt(gamma * pL / rL), abs(unR) + sqrt(gamma * pR / rR));

    for (int i = 0; i < 4; i++) {
        flux[i] = 0.5 * (FnL[i] + FnR[i] - lambdaMax * (UR[i] - UL[i]));
        flux[i] *= length;
    }
}

void Solver_2D::LLF_test(const double* UL, const double* UR, const double nx, const double ny, const double length, double* flux) {
    ////以下为原来的。目前存在发散的问题，可能是ULUR混淆
    //计算Fn、lambda
    const double& gamma = Constant::gamma;
    double Fn_L[4], Fn_R[4];//此处Fn为F·n
    double lambdaL, lambdaR;
    //计算Fn展开
    double ruvpL[4], ruvpR[4];
    Math_2D::U_2_ruvp(UL, ruvpL, gamma);
    Math_2D::U_2_ruvp(UR, ruvpR, gamma);
    //Math_2D::ruvp_2_Fn_lambda_2D(ruvp, Fn, lambda, nx, ny, gamma);
    double rhoL = ruvpL[0];
    double uL = ruvpL[1];
    double vL = ruvpL[2];
    double pL = ruvpL[3];
    double unL = uL * nx + vL * ny;
    double EL = Math_2D::get_E(ruvpL, gamma);
    Fn_L[0] = rhoL * unL;//rho*u*nx + rho*v*ny = rho*un
    Fn_L[1] = rhoL * uL * unL + pL * nx;//(rho*u^2+p)*nx + (rho*u*v)*ny = rho * u * un + p * nx
    Fn_L[2] = rhoL * vL * unL + pL + ny;//rho * v * un + p + ny
    Fn_L[3] = (rhoL * EL + pL) * unL;//(rho * E + p) * un
    lambdaL = abs(unL) + sqrt(gamma * pL / rhoL);

    double rhoR = ruvpR[0];
    double uR = ruvpR[1];
    double vR = ruvpR[2];
    double pR = ruvpR[3];
    double unR = uR * nx + vR * ny;
    double ER = Math_2D::get_E(ruvpR, gamma);
    Fn_R[0] = rhoR * unR;//rho*u*nx + rho*v*ny = rho*un
    Fn_R[1] = rhoR * uR * unR + pR * nx;//(rho*u^2+p)*nx + (rho*u*v)*ny = rho * u * un + p * nx
    Fn_R[2] = rhoR * vR * unR + pR + ny;//rho * v * un + p + ny
    Fn_R[3] = (rhoR * ER + pR) * unR;//(rho * E + p) * un
    lambdaR = abs(unR) + sqrt(gamma * pR / rhoR);

    //计算flux
    double lambdaMax = max(lambdaL, lambdaR);
    for (int i = 0; i < 4; i++) {
        flux[i] = 0.5 * (Fn_L[i] + Fn_R[i]) - 0.5 * lambdaMax * (UR[i] - UL[i]);
        flux[i] *= length;
    }


}

//void Solver_2D::calFlux_Roe_2() {
//    //功能：计算单元数值通量
//
//    //准备工作
//    FVM_2D* f = FVM_2D::pFVM2D;
//    for (int ie = 0; ie < f->elements.size(); ie++) {
//        for (int j = 0; j < 4; j++) {
//            f->elements[ie].Flux[j] = 0;//单元数值通量Flux清零，为后面加减做准备
//            if (global::flag_reconstruct == _REC_linear)f->elements[ie].updateSlope_Barth(f);//线性重构，即更新单元分布函数
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
//        //double NL = Neighbor(1, i);左右单元
//        //double NR = Neighbor(2, i);
//        //double DX = VECX(1, i);
//        //double DY = VECX(2, i);
//        //double n1 = FNode(1, i);节点坐标
//        //double n2 = FNode(2, i);
//
//
//        //double xmid, ymid;//边中点坐标
//        //pE->getxy(f, xmid, ymid);
//        //double sav1n, sav2n;//边法向？？？？？？？？？？？问问老师
//        ////pE->getDirectionN(sav1n, sav2n);
//        //double sav1 = (f->getNodeByID(pE->nodes[1])->x - f->getNodeByID(pE->nodes[0])->x);
//        //double sav2 = (f->getNodeByID(pE->nodes[1])->y - f->getNodeByID(pE->nodes[0])->y);
//        //double sav = sqrt(sav1 * sav1 + sav2 * sav2);
//        //sav1n = sav1 / sav;
//        //sav2n = sav2 / sav;
//
//        ////新增 直接导致发散
//        double sav1, sav2, sav1n, sav2n, sav;
//        pE->getDirectionN(sav1, sav2);
//        sav = sqrt(sav1 * sav1 + sav2 * sav2);
//        sav1n = sav1 / sav;
//        sav2n = sav2 / sav;
//
//
//        //do ig=1, NG 不知道干嘛的，难道是group？
//
//        //计算左ruvp
//        //!计算左单元在边中点处的ruvp
//        Eigen::Vector4d UL = pE->get_UL();
//        double _UL[4]{ UL[0],UL[1],UL[2],UL[3] };
//        double ruvp_L[4];
//        Math_2D::U_2_ruvp(_UL, ruvp_L, Constant::gamma);
//        double rl = ruvp_L[0];
//        double ul = ruvp_L[1];
//        double vl = ruvp_L[2];
//        double pl = ruvp_L[3];
//        //计算左单元中心点处的ruvp
//        double* UL_ele = pE->pElement_L->U;
//        double ruvp_L_ele[4];
//        Math_2D::U_2_ruvp(UL_ele, ruvp_L_ele, Constant::gamma);
//        double rlt = ruvp_L_ele[0];
//        double ult = ruvp_L_ele[1];
//        double vlt = ruvp_L_ele[2];
//        double plt = ruvp_L_ele[3];
//        //限制器。若边界处的值太离谱，则用单元中心处的值替代
//        if (rl < 0 || pl < 0 || rl>100000 || pl>1000000) {
//            rl = rlt;
//            ul = ult;
//            vl = vlt;
//            pl = plt;
//        }
//        double d_small = 1.0e-16;
//        rl = max(rl, d_small);
//        pl = max(pl, d_small);
//        rlt = max(rlt, d_small);
//        plt = max(plt, d_small);
//
//        //计算右ruvp
//        double rr;
//        double ur;
//        double vr;
//        double pr;
//        double rrt;
//        double urt;
//        double vrt;
//        double prt;
//        //判断边界类型
//        if (pE->pElement_R != nullptr) {//内部
//
//            //!计算右单元在边中点处的ruvp
//            Eigen::Vector4d UR = pE->get_UR();
//            double _UR[4]{ UR[0],UR[1],UR[2],UR[3] };
//            double ruvp_R[4];
//            Math_2D::U_2_ruvp(_UR, ruvp_R, Constant::gamma);
//            rr = ruvp_R[0];
//            ur = ruvp_R[1];
//            vr = ruvp_R[2];
//            pr = ruvp_R[3];
//            //计算右单元中心点处的ruvp
//            double* UR_ele = pE->pElement_R->U;
//            double ruvp_R_ele[4];
//            Math_2D::U_2_ruvp(UR_ele, ruvp_R_ele, Constant::gamma);
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
//            const int bType = f->boundaryManager.vBoundarySets[pE->setID - 1].type;
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
//                double x_edge, y_edge; pE->getxy(f, x_edge, y_edge);
//                std::vector<double> directionN = pE->getDirectionN(f);
//                double nx = directionN[0];
//                double ny = directionN[1];
//                Eigen::Matrix4d T = get_matrixT(nx, ny, 1);
//                Eigen::Matrix4d T_inv = get_matrixT(nx, ny, -1);
//                //计算U_L、U_R
//                Eigen::Vector4d U_L = pE->pElement_L->get_U(x_edge, y_edge);
//                Eigen::Vector4d U_R;
//                double* UR_ele;
//                {
//                    //找到这条edge对应的虚拟pElement_R(本来该edge没有pElement_R，但我们根据周期边界找到它的虚拟pElement_R)
//                    Edge_2D* pEdge_1 = f->boundaryManager.get_pairEdge_periodic(pE);//找到edge在周期边界中对应的edge_1
//                    Element_T3* virtual_pElement_R = pEdge_1->pElement_L;//pEdge_1的pElement_L即为pEdge的虚拟pElement_R
//                    UR_ele = virtual_pElement_R->U;
//                    //获取edge_1的中点坐标
//                    double x_virtual_edge, y_virtual_edge;
//                    pEdge_1->getxy(f, x_virtual_edge, y_virtual_edge);
//                    //重构。根据虚拟单元计算U_R
//                    U_R = virtual_pElement_R->get_U(x_virtual_edge, y_virtual_edge);
//                }
//
//                //计算右单元在边中点处的ruvp
//                Eigen::Vector4d UR = U_R;
//                double _UR[4]{ UR[0],UR[1],UR[2],UR[3] };
//                double ruvp_R[4];
//                Math_2D::U_2_ruvp(_UR, ruvp_R, Constant::gamma);
//                rr = ruvp_R[0];
//                ur = ruvp_R[1];
//                vr = ruvp_R[2];
//                pr = ruvp_R[3];
//                //计算右单元中心点处的ruvp
//                double ruvp_R_ele[4];
//                Math_2D::U_2_ruvp(UR_ele, ruvp_R_ele, Constant::gamma);
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
//        const double& gama = Constant::gamma;
//        const double ga1 = gama - 1;
//        double vaml = ult * ult + vlt * vlt;
//        double vamr = urt * urt + vrt * vrt;
//        double hlt = plt * gama / (rlt * (gama - 1.)) + 0.5 * vaml;
//        double hrt = prt * gama / (rrt * (gama - 1.)) + 0.5 * vamr;
//        double acl = sqrt(gama * plt / rlt);
//        double acr = sqrt(gama * prt / rrt);
//
//        ////Roe平均
//        double rrorl = sqrt(rrt / rlt);
//        double rrorlp1 = 1. + rrorl;
//        double rm = sqrt(rrt * rlt);
//        double um = (ult + urt * rrorl) / rrorlp1;
//        double vm = (vlt + vrt * rrorl) / rrorlp1;
//        double hm = (hlt + hrt * rrorl) / rrorlp1;
//        double vm2 = um * um + vm * vm;
//        double am2 = ga1 * abs(hm - 0.5 * vm2);
//        double am = sqrt(am2);
//        double am2i = 1. / am2;
//        double ami = 1. / am;
//        //特征值
//
//        double unormal = sav1 * um + sav2 * vm;
//        double anormal = am * sav;
//        double dacou = acr - acl;
//        double eig1 = abs(unormal);
//        double eig2 = abs(unormal + anormal);  //, deigen)
//        double eig3 = abs(unormal - anormal); //, deigen)
//
//
//        ////entropy fix, H
//        double yita = pE->pElement_L->deltaeig;//使用Compute_Deltaeig()进行初始化，该函数未完成
//        if (pE->pElement_R != nullptr) yita = max(yita, pE->pElement_R->deltaeig);
//        double epss = 2. * yita * sav;
//
//        //!!Harten entropy fix
//        if (eig1 < epss) eig1 = 0.5 * (eig1 * eig1 / epss + epss);
//        if (eig2 < epss) eig2 = 0.5 * (eig2 * eig2 / epss + epss);
//        if (eig3 < epss) eig3 = 0.5 * (eig3 * eig3 / epss + epss);
//
//        //计算hl hr
//        vaml = ul * ul + vl * vl;
//        vamr = ur * ur + vr * vr;
//        double hl = pl * gama / (rl * (gama - 1.)) + 0.5 * vaml;
//        double hr = pr * gama / (rr * (gama - 1.)) + 0.5 * vamr;
//
//        //inviscid flux
//        double dunormal = (sav1n * (ur - ul) + sav2n * (vr - vl));
//        double alph1 = eig1 * (rr - rl - (pr - pl) * am2i);
//        double alph2 = 0.5 * am2i * eig2 * (pr - pl + rm * am * dunormal);
//        double alph3 = 0.5 * am2i * eig3 * (pr - pl - rm * am * dunormal);
//        double alph4 = alph1 + alph2 + alph3;
//        double alph5 = am * (alph2 - alph3);
//        double alph6 = eig1 * rm * (ur - ul - sav1n * dunormal);
//        double alph7 = eig1 * rm * (vr - vl - sav2n * dunormal);
//
//        double ulnormal = sav1 * ul + sav2 * vl;
//        double urnormal = sav1 * ur + sav2 * vr;
//        double rulnormal = rl * ulnormal;
//        double rurnormal = rr * urnormal;
//        double unormaln = (unormal) / sav;
//
//        //数值通量
//        double flux[4];
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
//            std::string str = global::DoubleArray_2_String(f->edges[ie].pElement_L->U, 4) + "\n";
//            global::writeLogAndCout(str);
//        }
//    }
//
//}

void Solver_2D::Compute_Deltaeig() {
    //准备工作
    FVM_2D* f = FVM_2D::pFVM2D;
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].deltaeig = 0;//每一轮deltaeig清零
        }
    }
    //计算每个单元的deltaeig
    const double gamma = Constant::gamma;
    for (int ie = 0; ie < f->edges.size(); ie++) {
        Edge_2D* pE = &(f->edges[ie]);

        //double xmid, ymid;//边中点坐标
        //pE->getxy(f, xmid, ymid);
        //double sav1n, sav2n;//边法向？？？？？？？？？？？问问老师
        ////pE->getDirectionN(sav1n, sav2n);
        //double sav1 = (f->getNodeByID(pE->nodes[1])->x - f->getNodeByID(pE->nodes[0])->x);
        //double sav2 = (f->getNodeByID(pE->nodes[1])->y - f->getNodeByID(pE->nodes[0])->y);
        //double sav = sqrt(sav1 * sav1 + sav2 * sav2);
        //sav1n = sav1 / sav;
        //sav2n = sav2 / sav;

        ////新增 直接导致发散
        double sav1, sav2, sav1n, sav2n, sav;
        pE->getDirectionN(sav1, sav2);
        sav = sqrt(sav1 * sav1 + sav2 * sav2);
        sav1n = sav1 / sav;
        sav2n = sav2 / sav;


        //计算左ruvp
        //计算左单元在边中点处的ruvp
        double rl, ul, vl, pl, acl;
        double rr, ur, vr, pr, acr;
        if (pE->pElement_R != nullptr) {
            //计算左单元在边中点处的ruvp
            Eigen::Vector4d UL = pE->get_UL();
            double _UL[4]{ UL[0],UL[1],UL[2],UL[3] };
            double ruvp_L[4];
            Math_2D::U_2_ruvp(_UL, ruvp_L, Constant::gamma);
            rl = ruvp_L[0];
            ul = ruvp_L[1];
            vl = ruvp_L[2];
            pl = ruvp_L[3];
            acl = sqrt(gamma * pl / rl);
            //计算右单元在边中点处的ruvp
            Eigen::Vector4d UR = pE->get_UR();
            double _UR[4]{ UR[0],UR[1],UR[2],UR[3] };
            double ruvp_R[4];
            Math_2D::U_2_ruvp(_UR, ruvp_R, Constant::gamma);
            rr = ruvp_R[0];
            ur = ruvp_R[1];
            vr = ruvp_R[2];
            pr = ruvp_R[3];
            acr = sqrt(gamma * pr / rr);
            //计算deltaeig
            double deig = abs(sav1n * (ul - ur) + sav2n * (vl - vr)) + abs(acl - acr);
            pE->pElement_L->deltaeig = max(pE->pElement_L->deltaeig, deig);
            pE->pElement_R->deltaeig = max(pE->pElement_R->deltaeig, deig);
        }
    }
}

void Solver_2D::cal_ruvp_farfield_new(const double nx, const double ny, double* ruvp, const double* ruvp_inf) {
    //边界元的Element_R为null，因此nxny一定朝外
    const double rho = ruvp[0];
    const double u = ruvp[1];
    const double v = ruvp[2];
    const double p = ruvp[3];
    const double rho_inf = ruvp_inf[0];
    const double u_inf = ruvp_inf[1];
    const double v_inf = ruvp_inf[2];
    const double p_inf = ruvp_inf[3];
    const double gamma = Constant::gamma;

    //坐标变换
    const double rho_n = rho;
    const double u_n = u * nx + v * ny;//u cost + v sint
    const double v_n = -u * ny + v * nx;//此处有修改。法线向外为正 - u-sint + v cost
    //const double v_n = u * ny - v * nx;
    const double p_n = p;
    const double rho_n_inf = rho_inf;
    const double p_n_inf = p_inf;
    const double u_n_inf = u_inf * nx + v_inf * ny;
    const double v_n_inf = -u_inf * ny + v_inf * nx;//此处有修改

    //if (rho_n < 0)system("pause");
    double a_n2 = gamma * p_n / rho_n;
    if (a_n2 < 0) { std::cout << "Error: a_n2 < 0\n";  system("pause"); }
    const double a_n = sqrt(a_n2);
    double a_n_inf2 = gamma * p_n_inf / rho_n_inf;
    if (a_n_inf2 < 0) { std::cout << "Error: a_n_inf2 < 0\n";  system("pause"); }
    const double a_n_inf = sqrt(a_n_inf2);

    double v_n_tmp;
    double u_n_tmp;
    double rho_n_tmp;
    double p_n_tmp;

    if (u_n_inf * u_n_inf < a_n_inf * a_n_inf) {//|Vn_inf|<|a_inf|，亚音速
        double tmp_a = (gamma - 1.0) / 4.0 * (u_n - u_n_inf) + 0.5 * (a_n + a_n_inf);// 1/2 *( (ga-1)/2 * (u-uinf) + (a+ainf) )
        double tmp_a2 = tmp_a * tmp_a;
        if (u_n_inf <= 0.0) {//亚音速入口，提3个边界条件
            double tmp_s_inf = p_n_inf / pow(rho_n_inf, gamma);//
            rho_n_tmp = pow(tmp_a2 / tmp_s_inf / gamma, 1.0 / (gamma - 1.0));//边
            u_n_tmp = 0.5 * (u_n + u_n_inf + 2.0 * a_n / (gamma - 1.0) - 2.0 * a_n_inf / (gamma - 1.0));
            v_n_tmp = v_n_inf;//边
            p_n_tmp = tmp_a2 * rho_n_tmp / gamma;//边
        }
        else {//亚音速出口，提1个边界条件 u_n_tmp
            double tmp_s = p_n / pow(rho_n, gamma);//内 //pow(x,y):若x<0且y非整数，或者x=0且y<=0，将出现结果错误
            rho_n_tmp = pow(tmp_a2 / tmp_s / gamma, 1.0 / (gamma - 1.0));//内
            u_n_tmp = 0.5 * (u_n + u_n_inf + 2.0 * a_n / (gamma - 1.0) - 2.0 * a_n_inf / (gamma - 1.0));//边 //法向按公式计算u_n_i_m_η后，与u_n平均
            v_n_tmp = v_n;
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
