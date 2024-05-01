#include "CIteration.h"

#include "Evolve.h"
#include "../space/restrict/Restrict.h"
#include "../global/GlobalPara.h"
#include "CalculateDt.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"
#include "../output/LogWriter.h"
#include "EvolveGPU.h"


void CIteration::iteration_host_2(myfloat& t, myfloat T, GPU::GPUSolver2* s) {
    update_host_ruvp_Uold(s->elementField_host.U, s->element_vruvp, s->element_U_old, s->elementField_host.num, GlobalPara::constant::gamma);
    m_EvolveHost_2_addRK3(t, T, s->element_host, s->elementField_host, s->edge_host, s->node_host, s->element_vruvp);
    U2NITS::Space::Restrict::modifyElementFieldU2d(s->element_host, s->elementField_host);

}

void CIteration::iteration_useGPU(myfloat& t, myfloat T, GPU::GPUSolver2* s) {
    update_host_ruvp_Uold(s->elementField_host.U, s->element_vruvp, s->element_U_old, s->elementField_host.num, GlobalPara::constant::gamma);

    bool is_steady = GlobalPara::time::is_steady;
    bool is_explicit = GlobalPara::time::is_explicit;
    int temporal_choice = GlobalPara::time::time_advance;
    auto func_unimplemented = [](std::string str = "")->void {
        LogWriter::logAndPrintError("unimplemented " + str + " @CIteration::m_EvolveGPU_2\n");
        exit(-1);
    };

    // ����
    if (is_steady) {
        func_unimplemented("steady");
    }
    // �Ƕ���
    else {
        // ��������ʱ�䲽��
        myfloat dt = U2NITS::Time::calculateGlobalDt(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::time::CFL,
            GlobalPara::constant::R, s->element_host, s->edge_host, s->element_vruvp);
        // �޸ĵ�ǰ����ʱ��t
        CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
        pCDS->set_t(t);
        pCDS->set_dt(0.0);
        t += dt;

        GPU::ElementFieldSoA::cuda_memcpy(&s->elementField_device, &s->elementField_host, cudaMemcpyHostToDevice);

        switch (temporal_choice) {
        case _EVO_explicit:
            GPU::Time::evolveSingleStep_device(dt, s->element_device, s->node_device, s->edge_device, s->elementField_device, s->boundary_device, s->sDevicePara);
            break;
        case _EVO_rk3:
            func_unimplemented("rk3");
            break;
        default:
            func_unimplemented("temporal_choice");
            break;
        }
    }


}

void CIteration::update_host_ruvp_Uold(myfloat* U[4], myfloat* ruvp[4], myfloat* U_old[4], myint num, myfloat gamma) {
    /*
    ��U����element_vruvp.element_vruvp���ڼ���ʱ�䲽��Dt�Լ��������
    U:rho,rho_u,rho_v,rho_E ruvp:rho,u,v,p
    */
    for (int i = 0; i < num; i++) {
        // ����host�˵�U_old
        for (int j = 0; j < 4; j++) {
            U_old[j][i] = U[j][i];
        }
        // ����host�˵�ruvp  
        ruvp[0][i] = U[0][i];
        ruvp[1][i] = U[1][i] / U[0][i];
        ruvp[2][i] = U[2][i] / U[0][i];
        ruvp[3][i] = ruvp[0][i] * (gamma - 1) * (U[3][i] / U[0][i] - (ruvp[1][i] * ruvp[1][i] + ruvp[2][i] * ruvp[2][i]) * 0.5);
    }

}

void CIteration::m_EvolveHost_2_addRK3(myfloat& physicalTime, myfloat maxPhysicalTime, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, GPU::NodeSoA& node_host, myfloat* element_vruvp[4]) {
    /*
    ʱ���ƽ�
    ���������þֲ�ʱ���ƽ���t��T��Ч
    �Ƕ�������������ʱ���ƽ�����Ҫ�ı�t��ֵ�����ȼ�����������dt��Ȼ���ƽ����ɲ���RK1��RK3
    */
    using namespace U2NITS::Time;
    bool is_steady = GlobalPara::time::is_steady;
    bool is_explicit = GlobalPara::time::is_explicit;
    int temporal_choice = GlobalPara::time::time_advance;

    // ����
    if (is_steady) {
        physicalTime += 1e-3;// �ı�t��Ϊ�����tecplotʱ�ܹ������еķ�ʽ��ȡ
        evolveSteadyLocalTimeStep(element_host, node_host, edge_host, elementField_host, element_vruvp);
    }
    // �Ƕ���
    else {
        // ��������ʱ�䲽��
        myfloat dt = U2NITS::Time::calculateGlobalDt(physicalTime, maxPhysicalTime, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::time::CFL,
            GlobalPara::constant::R, element_host, edge_host, element_vruvp);
        // �޸ĵ�ǰ����ʱ��t
        CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
        pCDS->set_t(physicalTime);
        pCDS->set_dt(0.0);
        physicalTime += dt;

        switch (temporal_choice) {
        case _EVO_explicit:
            evolveSingleStep(dt, element_host, node_host, edge_host, elementField_host);
            break;
        case _EVO_rk3:
            evolveRungeKutta3(dt, element_host, node_host, edge_host, elementField_host);
            break;
        default:
            LogWriter::logAndPrintError("invalid temporal_choice @U2NITS::Time::EvolveHost_2.\n");
            exit(-1);
        }
    }
}

