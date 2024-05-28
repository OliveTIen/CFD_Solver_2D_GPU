#include "Solver_2D.h"
#include "FVM_2D.h"
#include "output/LogWriter.h"
#include "space/RiemannSolver.h"
#include "space/Reconstructor.h"

//myfloat Solver_2D::RK3alpha[6]{ 0.0,1.0 / 4.0,1.0 / 6.0,3.0 / 8.0 ,0.5,1.0 };
//myfloat Solver_2D::RK5alpha[6]{ 0.0,1.0 / 4.0,1.0 / 6.0,3.0 / 8.0 ,0.5,1.0 };


//Eigen::Matrix4d Solver_2D::get_matrixT(myfloat nx, myfloat ny, int flag) {
//    Eigen::Matrix4d matT;
//    ny *= flag;//��������sin��Ϊ�෴��
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
    //��ֵͨ��
    calFlux();
    //ʱ���ƽ�
    FVM_2D* f = FVM_2D::getInstance();
    myfloat omega;//���
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
    //����(�ֲ�)�ṹ�壬���ڴ洢Uֵ
    struct Data {
        myfloat U[4]{};
    };
    std::vector<Data> elements_old(f->elements.size());
    std::vector<Data> k1(elements_old), k2(k1), k3(k1);
    //����
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            elements_old[ie].U[j] = f->elements[ie].U[j];
        }
    }
    
    //����������洢����ֹ��������
    std::vector<myfloat>areas; areas.resize(f->elements.size());
    for (int ie = 0; ie < f->elements.size(); ie++) {
        areas[ie] = f->elements[ie].calArea(f);
    }

    
    //1.k1=f(t_n,y_n)
    //1.1.����Uֵ������ֵͨ��
    calFlux();
    //1.2.��k1
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            k1[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / areas[ie];//ע�⸺��
        }
    }
    //2.k2=f(t_n+dt/2,y_n+dt/2*k1)
    //2.1.��Uֵy_n+dt/2*k1
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].U[j] = elements_old[ie].U[j] + dt/2.0 * k1[ie].U[j];
        }
    }
    //2.2.����Uֵ������ֵͨ��
    calFlux();
    //2.3.��k2
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            k2[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / areas[ie];//ע�⸺��
        }
    }
    //3.k3=f(t_n+dt,y_n-dt*k1+2*dt*k2)
    //3.1.��Uֵy_n-dt*k1+2*dt*k2
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].U[j] = elements_old[ie].U[j] - dt * k1[ie].U[j] + 2.0 * dt * k2[ie].U[j];
        }
    }
    //3.2.����Uֵ������ֵͨ��
    calFlux();
    //3.3.��k3
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            k3[ie].U[j] = -1.0 * f->elements[ie].Flux[j] / areas[ie];//ע�⸺��
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
    �ѳ�ʼ���Ĳ�����
        ��Ԫ����

    ����ͨ�����̣�
	1. ��ʼ�����е�Ԫ��ֵͨ��Ϊ0
    2. ���㵥Ԫ�ݶ�
        ���ݵ�Ԫ���ꡢ��ԪU���ھ����ꡢ�ھ�U�������ݶ�
        �����������ݶ�
        ����쳣ֵ
        ���룺��Ԫ���ꡢ��ԪU���ھ����ꡢ�ھ�U
        ������ݶ�
    3. ����߽���ֵͨ��
        ����߷�������������
        ��������ҵ�Ԫ�ڱ߽��Uֵ(���Ҽ���UL��UR)
        ���ݱ����ͣ�����UL��UR
        ������������
            ���룺UL UR nx ny edgeLength
            ������߽����ֵͨ��
    4. ���㵥Ԫ��ֵͨ��
        ���߽���ֵͨ���ӵ���Ԫ��ֵͨ����
    5. ��ʽʱ���ƽ�

    */


    FVM_2D* f = FVM_2D::getInstance();

//#pragma omp parallel for
    //׼����������ʼ����Ԫ��ֵͨ�����ֲ�����
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].Flux[j] = 0;//��Ԫ��ֵͨ��Flux���㣬Ϊ����Ӽ���׼��
            f->elements[ie].deltaeig = 0;//ÿһ��deltaeig����
            if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_constant
                && GlobalPara::physicsModel::equation == _EQ_euler) {
                // ����ǳ����ع�+Euler������������ݶ�
            }
            else {
                // ������Ҫ�����ݶȣ���������������
                Reconstructor::Element_T3_updateSlope_Barth(f, &(f->elements[ie]));
                //f->elements[ie].updateSlope_Barth(f);//�����ع���Barth������
            }
        }
    }

//#pragma omp parallel for
    //ÿ���߼�����ճͨ����Ȼ����ݷ���ֱ�Ӽ������൥Ԫ��Flux�����б߱��������е�Ԫ��FluxҲ�ͼ��������
    for (int iedge = 0; iedge < f->edges.size(); iedge++) {
        ////����ÿ���ߵ�����ͨ��
        myfloat flux[4]{};

        {
            FVM_2D* f = FVM_2D::getInstance();
            Edge_2D* pE = &(f->edges[iedge]);
            if (pE->pElement_R == nullptr) {//�߽�
                const int bType = f->boundaryManager.boundaries[pE->setID - 1].type;
                switch (bType) {
                    //�Գơ���ŷ�����̶��ԣ��ԳƱ߽缴��ճ�̱ڱ߽�
                case _BC_symmetry:
                    getEdgeFlux_wallNonViscous(pE, flux);
                    break;

                    //��ճ�̱�
                case _BC_wall_nonViscous:
                    getEdgeFlux_wallNonViscous(pE, flux);
                    break;
                    //�޻��ƾ���
                case _BC_wall_adiabat:
                    
                    break;
                    //���
                case _BC_inlet:
                    getEdgeFlux_farfield(pE, GlobalPara::boundaryCondition::_2D::inlet::ruvp, flux);
                    break;
                    //����
                case _BC_outlet:
                    getEdgeFlux_farfield(pE, GlobalPara::boundaryCondition::_2D::outlet::ruvp, flux);
                    break;
                    //Զ��
                case _BC_inf:
                    getEdgeFlux_farfield(pE, GlobalPara::boundaryCondition::_2D::inf::ruvp, flux);
                    break;
                }
                //����
                if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) getEdgeFlux_periodic(pE, flux);//���10�����ڱ߽�
            }
            else//�ڲ�
            {
                getEdgeFlux_inner(pE, flux);
            }

        }

        ////�������൥Ԫ��Flux
        for (int j = 0; j < 4; j++) {
            if (f->edges[iedge].pElement_R != nullptr) {//�ڲ�edge
                f->edges[iedge].pElement_L->Flux[j] += flux[j];
                f->edges[iedge].pElement_R->Flux[j] -= flux[j];
            }
            else {//�߽�edge
                f->edges[iedge].pElement_L->Flux[j] += flux[j];
            }
        }
    }

}

void Solver_2D::getEdgeFlux_inner(Edge_2D* pE, myfloat* flux) {
    //���ܣ������ڲ��߽����ֵͨ��
    
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

    //����U_R
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
    //�������ܣ�����Զ���߽���������߽���ֵͨ��
    FVM_2D* f = FVM_2D::getInstance();
    myfloat nx, ny; 
    myfloat x_edge, y_edge;
    myfloat U_L[4];
    myfloat ruvp_L[4];
    //����ruvp_L
    pE->getDirectionN(nx, ny);
    pE->getxy(FVM_2D::getInstance(), x_edge, y_edge);
    pE->pElement_L->get_U(x_edge, y_edge, U_L);
    Math_2D::U_2_ruvp(U_L, ruvp_L, GlobalPara::constant::gamma);
    //���ruvp_L�Ϸ��� ����Ƿ�ɢ
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
    //����Զ���߽���������ruvp_L
    modify_ruvpL_in_farfield(nx, ny, ruvp_L, ruvp_inf);
    Math_2D::ruvp_2_U(ruvp_L, U_L, GlobalPara::constant::gamma);
    RiemannSolve(U_L, U_L, nx, ny, pE->getLength(), flux,
        GlobalPara::inviscid_flux_method::flux_conservation_scheme);

}

void Solver_2D::getEdgeFlux_periodic(Edge_2D* pE, myfloat* flux) {
    //FVM_2D* f = FVM_2D::getInstance();
    ////������е����ꡢ����任����������任���󡢱߷��߷���
    //myfloat x_edge, y_edge; pE->getxy(f, x_edge, y_edge);
    //myfloat nx, ny;    pE->getDirectionN(nx, ny);
    ////����U_L��U_R
    //Eigen::Vector4d U_L = pE->pElement_L->get_U(x_edge, y_edge);
    //Eigen::Vector4d U_R;
    //{
    //    //�ҵ�����edge��Ӧ������pElement_R(������edgeû��pElement_R�������Ǹ������ڱ߽��ҵ���������pElement_R)
    //    Edge_2D* pE_pair = f->boundaryManager.get_pairEdge_periodic(pE);//�ҵ�edge�����ڱ߽��ж�Ӧ��edge_1
    //    Element_2D* pElement_R = pE_pair->pElement_L;//pEdge_1��pElement_L��ΪpEdge������pElement_R
    //    //��ȡedge_1���е�����
    //    myfloat x_virtual_edge, y_virtual_edge;
    //    pE_pair->getxy(f, x_virtual_edge, y_virtual_edge);
    //    //�ع����������ⵥԪ����U_R
    //    U_R = pElement_R->get_U(x_virtual_edge, y_virtual_edge);
    //}
    ////������ֵͨ����������Ԫ��Ϊ���غ�������ֵͨ��
    //myfloat UL[4]{ U_L[0],U_L[1],U_L[2],U_L[3] };
    //myfloat UR[4]{ U_R[0],U_R[1],U_R[2],U_R[3] };
    //RiemannSolve(UL, UR, nx, ny, pE->getLength(), flux,
    //    GlobalPara::inviscid_flux_method::flux_conservation_scheme);

    // ����д
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
//    //���ܣ����㵥Ԫ��ֵͨ��
//
//    //׼������
//    FVM_2D* f = FVM_2D::getInstance();
//    for (int ie = 0; ie < f->elements.size(); ie++) {
//        for (int j = 0; j < 4; j++) {
//            f->elements[ie].Flux[j] = 0;//��Ԫ��ֵͨ��Flux���㣬Ϊ����Ӽ���׼��
//            if (GlobalStatic::flag_reconstruct == _REC_linear)f->elements[ie].updateSlope_Barth(f);//�����ع��������µ�Ԫ�ֲ�����
//        }
//    }
//
//    //��Ԫdeltaeig���㣬Ȼ��ֵ
//    Compute_Deltaeig();
//
//    for (int ie = 0; ie < f->edges.size(); ie++) {
//        Edge_2D* pE = &(f->edges[ie]);
//
//        //2023-8-16 �·���
//        //myfloat NL = Neighbor(1, i);���ҵ�Ԫ
//        //myfloat NR = Neighbor(2, i);
//        //myfloat DX = VECX(1, i);
//        //myfloat DY = VECX(2, i);
//        //myfloat n1 = FNode(1, i);�ڵ�����
//        //myfloat n2 = FNode(2, i);
//
//
//        //myfloat xmid, ymid;//���е�����
//        //pE->getxy(f, xmid, ymid);
//        //myfloat sav1n, sav2n;//�߷��򣿣�������������������������ʦ
//        ////pE->getDirectionN(sav1n, sav2n);
//        //myfloat sav1 = (f->getNodeByID(pE->nodes[1])->x - f->getNodeByID(pE->nodes[0])->x);
//        //myfloat sav2 = (f->getNodeByID(pE->nodes[1])->y - f->getNodeByID(pE->nodes[0])->y);
//        //myfloat sav = sqrt(sav1 * sav1 + sav2 * sav2);
//        //sav1n = sav1 / sav;
//        //sav2n = sav2 / sav;
//
//        ////���� ֱ�ӵ��·�ɢ
//        myfloat sav1, sav2, sav1n, sav2n, sav;
//        pE->getDirectionN(sav1, sav2);
//        sav = sqrt(sav1 * sav1 + sav2 * sav2);
//        sav1n = sav1 / sav;
//        sav2n = sav2 / sav;
//
//
//        //do ig=1, NG ��˹���ֵ�
//
//        //������ruvp
//        //!������Ԫ�ڱ��е㴦��ruvp
//        Eigen::Vector4d UL = pE->get_UL();
//        myfloat _UL[4]{ UL[0],UL[1],UL[2],UL[3] };
//        myfloat ruvp_L[4];
//        Math_2D::U_2_ruvp(_UL, ruvp_L, GlobalPara::constant::gamma);
//        myfloat rl = ruvp_L[0];
//        myfloat ul = ruvp_L[1];
//        myfloat vl = ruvp_L[2];
//        myfloat pl = ruvp_L[3];
//        //������Ԫ���ĵ㴦��ruvp
//        myfloat* UL_ele = pE->pElement_L->U;
//        myfloat ruvp_L_ele[4];
//        Math_2D::U_2_ruvp(UL_ele, ruvp_L_ele, GlobalPara::constant::gamma);
//        myfloat rlt = ruvp_L_ele[0];
//        myfloat ult = ruvp_L_ele[1];
//        myfloat vlt = ruvp_L_ele[2];
//        myfloat plt = ruvp_L_ele[3];
//        //�����������߽紦��ֵ̫���ף����õ�Ԫ���Ĵ���ֵ���
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
//        //������ruvp
//        myfloat rr;
//        myfloat ur;
//        myfloat vr;
//        myfloat pr;
//        myfloat rrt;
//        myfloat urt;
//        myfloat vrt;
//        myfloat prt;
//        //�жϱ߽�����
//        if (pE->pElement_R != nullptr) {//�ڲ�
//
//            //!�����ҵ�Ԫ�ڱ��е㴦��ruvp
//            Eigen::Vector4d UR = pE->get_UR();
//            myfloat _UR[4]{ UR[0],UR[1],UR[2],UR[3] };
//            myfloat ruvp_R[4];
//            Math_2D::U_2_ruvp(_UR, ruvp_R, GlobalPara::constant::gamma);
//            rr = ruvp_R[0];
//            ur = ruvp_R[1];
//            vr = ruvp_R[2];
//            pr = ruvp_R[3];
//            //�����ҵ�Ԫ���ĵ㴦��ruvp
//            myfloat* UR_ele = pE->pElement_R->U;
//            myfloat ruvp_R_ele[4];
//            Math_2D::U_2_ruvp(UR_ele, ruvp_R_ele, GlobalPara::constant::gamma);
//            rrt = ruvp_R_ele[0];
//            urt = ruvp_R_ele[1];
//            vrt = ruvp_R_ele[2];
//            prt = ruvp_R_ele[3];
//            //������
//            if (rr < 0 || pr < 0 || rr>100000 || pr>1000000) {
//                rr = rrt;
//                ur = urt;
//                vr = vrt;
//                pr = prt;
//            }
//        }
//
//        else {//�߽� ��δ���
//            const int bType = f->boundaryManager.boundaries[pE->setID - 1].type;
//            switch (bType) {
//                //�Գơ���ŷ�����̶��ԣ��ԳƱ߽缴��ճ�̱ڱ߽�
//            case _BC_symmetry:
//                //����break
//
//                //��ճ�̱�
//            case _BC_wall_nonViscous:
//            {
//
//            }
//            //�޻��ƾ���
//            case _BC_wall_adiabat:
//                
//                //���
//            case _BC_inlet:
//                //return calEdgeFlux_LLF_farfield(pE, GlobalPara::boundaryCondition::_2D::inlet::ruvp);
//                //����
//            case _BC_outlet:
//                //return calEdgeFlux_LLF_farfield(pE, GlobalPara::boundaryCondition::_2D::outlet::ruvp);
//                //Զ��
//            case _BC_inf:
//                //return calEdgeFlux_LLF_farfield(pE, GlobalPara::boundaryCondition::_2D::inf::ruvp);
//                std::cout << "Error: δ��ɴ˲��ִ���(Solver_2D::calFlux_Roe_2)\n";
//                break;
//
//            default:
//                break;
//            }
//            //����
//            if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {//���10�����ڱ߽�
//                //������е����ꡢ����任����������任���󡢱߷��߷���
//                myfloat x_edge, y_edge; pE->getxy(f, x_edge, y_edge);
//                std::vector<myfloat> directionN = pE->getDirectionN(f);
//                myfloat nx = directionN[0];
//                myfloat ny = directionN[1];
//                Eigen::Matrix4d T = get_matrixT(nx, ny, 1);
//                Eigen::Matrix4d T_inv = get_matrixT(nx, ny, -1);
//                //����U_L��U_R
//                Eigen::Vector4d U_L = pE->pElement_L->get_U(x_edge, y_edge);
//                Eigen::Vector4d U_R;
//                myfloat* UR_ele;
//                {
//                    //�ҵ�����edge��Ӧ������pElement_R(������edgeû��pElement_R�������Ǹ������ڱ߽��ҵ���������pElement_R)
//                    Edge_2D* pEdge_1 = f->boundaryManager.get_pairEdge_periodic(pE);//�ҵ�edge�����ڱ߽��ж�Ӧ��edge_1
//                    Element_2D* virtual_pElement_R = pEdge_1->pElement_L;//pEdge_1��pElement_L��ΪpEdge������pElement_R
//                    UR_ele = virtual_pElement_R->U;
//                    //��ȡedge_1���е�����
//                    myfloat x_virtual_edge, y_virtual_edge;
//                    pEdge_1->getxy(f, x_virtual_edge, y_virtual_edge);
//                    //�ع����������ⵥԪ����U_R
//                    U_R = virtual_pElement_R->get_U(x_virtual_edge, y_virtual_edge);
//                }
//
//                //�����ҵ�Ԫ�ڱ��е㴦��ruvp
//                Eigen::Vector4d UR = U_R;
//                myfloat _UR[4]{ UR[0],UR[1],UR[2],UR[3] };
//                myfloat ruvp_R[4];
//                Math_2D::U_2_ruvp(_UR, ruvp_R, GlobalPara::constant::gamma);
//                rr = ruvp_R[0];
//                ur = ruvp_R[1];
//                vr = ruvp_R[2];
//                pr = ruvp_R[3];
//                //�����ҵ�Ԫ���ĵ㴦��ruvp
//                myfloat ruvp_R_ele[4];
//                Math_2D::U_2_ruvp(UR_ele, ruvp_R_ele, GlobalPara::constant::gamma);
//                rrt = ruvp_R_ele[0];
//                urt = ruvp_R_ele[1];
//                vrt = ruvp_R_ele[2];
//                prt = ruvp_R_ele[3];
//                //������
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
//        ////Roeƽ��
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
//        //����ֵ
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
//        myfloat yita = pE->pElement_L->deltaeig;//ʹ��Compute_Deltaeig()���г�ʼ�����ú���δ���
//        if (pE->pElement_R != nullptr) yita = max(yita, pE->pElement_R->deltaeig);
//        myfloat epss = 2. * yita * sav;
//
//        //!!Harten entropy fix
//        if (eig1 < epss) eig1 = 0.5 * (eig1 * eig1 / epss + epss);
//        if (eig2 < epss) eig2 = 0.5 * (eig2 * eig2 / epss + epss);
//        if (eig3 < epss) eig3 = 0.5 * (eig3 * eig3 / epss + epss);
//
//        //����hl hr
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
//        //��ֵͨ��
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
//        //���ߵ���ֵͨ���Ӽ�����Ԫ��
//        for (int j = 0; j < 4; j++) {
//            pE->pElement_L->Flux[j] += flux[j];
//            if (pE->pElement_R != nullptr) {//�ڲ�edge
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
    //׼������
    FVM_2D* f = FVM_2D::getInstance();
    for (int ie = 0; ie < f->elements.size(); ie++) {
        for (int j = 0; j < 4; j++) {
            f->elements[ie].deltaeig = 0;//ÿһ��deltaeig����
        }
    }
    //����ÿ����Ԫ��deltaeig
    const myfloat gamma = GlobalPara::constant::gamma;
    for (int ie = 0; ie < f->edges.size(); ie++) {
        Edge_2D* pE = &(f->edges[ie]);

        //myfloat xmid, ymid;//���е�����
        //pE->getxy(f, xmid, ymid);
        //myfloat sav1n, sav2n;//�߷��򣿣�������������������������ʦ
        ////pE->getDirectionN(sav1n, sav2n);
        //myfloat sav1 = (f->getNodeByID(pE->nodes[1])->x - f->getNodeByID(pE->nodes[0])->x);
        //myfloat sav2 = (f->getNodeByID(pE->nodes[1])->y - f->getNodeByID(pE->nodes[0])->y);
        //myfloat sav = sqrt(sav1 * sav1 + sav2 * sav2);
        //sav1n = sav1 / sav;
        //sav2n = sav2 / sav;

        ////���� ֱ�ӵ��·�ɢ
        myfloat sav1, sav2, sav1n, sav2n, sav;
        pE->getDirectionN(sav1, sav2);
        sav = sqrt(sav1 * sav1 + sav2 * sav2);
        sav1n = sav1 / sav;
        sav2n = sav2 / sav;


        //������ruvp
        //������Ԫ�ڱ��е㴦��ruvp
        myfloat rl, ul, vl, pl, acl;
        myfloat rr, ur, vr, pr, acr;
        if (pE->pElement_R != nullptr) {
            //���㵥Ԫ�ڱ��е㴦��ruvp
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
            //����deltaeig
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
    //�߽�Ԫ��Element_RΪnull�����nxnyһ������
    const myfloat rho = ruvp[0];
    const myfloat u = ruvp[1];
    const myfloat v = ruvp[2];
    const myfloat p = ruvp[3];
    const myfloat rho_inf = ruvp_inf[0];
    const myfloat u_inf = ruvp_inf[1];
    const myfloat v_inf = ruvp_inf[2];
    const myfloat p_inf = ruvp_inf[3];
    const myfloat gamma = GlobalPara::constant::gamma;

    //����任
    const myfloat rho_n = rho;
    const myfloat u_n = u * nx + v * ny;//u cost + v sint
    const myfloat v_n = -u * ny + v * nx;//�˴����޸ġ���������Ϊ�� - u sint + v cost
    //const myfloat v_n = u * ny - v * nx;
    const myfloat p_n = p;
    const myfloat rho_n_inf = rho_inf;
    const myfloat p_n_inf = p_inf;
    const myfloat u_n_inf = u_inf * nx + v_inf * ny;
    const myfloat v_n_inf = -u_inf * ny + v_inf * nx;//�˴����޸�

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

void Solver_2D::flux_viscous(Edge_2D* pE, myfloat* flux) {
    // Ӧ�Ž�Solver_2D::calFlux()��for each face boundary��
    // �˴�����ͬ�߽������Ĵ���ȫ�Ž�һ��������
    // edge�ֳ�Ϊface



    FVM_2D* pFVM2D = FVM_2D::getInstance();
    myfloat nx, ny;
    pE->getDirectionN(nx, ny);
    myfloat x_edge, y_edge;
    pE->getxy(pFVM2D, x_edge, y_edge);
    const myfloat& gamma = GlobalPara::constant::gamma;// ���ã�����ռ�ռ�
    const int nVar = 4;
    myfloat sav = pE->getRefLength();// ! �˴��Ѿ��޸ģ�ǰ����ݶ���������ҲҪ�޸ģ�
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
    myfloat rl{};// ��ʼ��Ϊ0.0
    myfloat ul{};// myfloat rl(ruvp_L[0])�൱�� myfloat rl = ruvp_L[0]
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
    myfloat rr{};// ��ʼ��Ϊ0.0
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

        // �̱�
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
            // !!!! ע��WangQ��"flux_viscous.f90��142��"�����õ�ruvp���õ�puvt��Ҫ��������Ϊ�����ɢ��t�ĺ�ɢ
            tgd[0][k] = 0.5 * (tgd[0][k] + pE->pElement_R->Ux[k] * dx)
                + yita * sav1n * (ruvp_R[k] - ruvp_L[k]) / dxk;
            tgd[1][k] = 0.5 * (tgd[1][k] + pE->pElement_R->Uy[k] * dy)
                + yita * sav2n * (ruvp_R[k] - ruvp_L[k]) / dxk;
        }
    }
    // ��������WangQ��"flux_viscous.f90��150��
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
    // reflen �����ٻ��ĳ��� ȡ���Ǿ���ֵ
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

        // �̱�
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
    // ����֤��air��Cp Cv��������������ϵʽ

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
