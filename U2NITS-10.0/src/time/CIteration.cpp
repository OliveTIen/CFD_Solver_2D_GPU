#include "CIteration.h"

#include "Evolve.h"
#include "../space/restrict/Restrict.h"
#include "../global/GlobalPara.h"
#include "CalculateDt.h"
#include "CalculateDtGPU.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"
#include "../output/LogWriter.h"
#include "EvolveGPU.h"
#include "../math/PhysicalKernelGPU.h"


void constPointerTest() {
    /*
    ѧϰconst��ָ�����η���ʹ�á����README(05-17)
    */
    const int a = 10;
    int b = 20;
    int* p = (int*) & a;

    const int* const q = &a;
    int* const r = &b;
    *r = 1;
}

void CIteration::iteration_host(myfloat& t, myfloat T, GPU::GPUSolver2* s) {

    update_host_ruvp_Uold(s->elementField_host.U, s->elementField_host.ruvp, s->elementField_host.Uold, s->elementField_host.num, GlobalPara::constant::gamma);

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
        t += 1e-3;// �ı�t��Ϊ�����tecplotʱ�ܹ������еķ�ʽ��ȡ������۲춯���Խ��е���
        evolve_steady_explicit_localTimeStep(s->element_host, s->node_host, s->edge_host, s->elementField_host, s->elementField_host.ruvp);
    }
    // �Ƕ���
    else {
        // ��������ʱ�䲽��
        myfloat dt = U2NITS::Time::get_global_dt_host(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Pr, GlobalPara::time::CFL,
            GlobalPara::constant::R, s->elementFieldVariable_dt_host, s->element_host, s->edge_host, s->elementField_host.ruvp, GlobalPara::physicsModel::equation);
        // ˫��շ����޸ı߽�
        CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
        pCDS->set_t(t);
        pCDS->set_dt(0.0);
        //
        t += dt;

        switch (temporal_choice) {
        case _EVO_explicit:
            evolve_unsteady_explicit(dt, s->element_host, s->node_host, s->edge_host, s->elementField_host);
            break;
        case _EVO_rk3:
            evolve_unsteady_rk3(dt, s->element_host, s->node_host, s->edge_host, s->elementField_host);
            break;
        default:
            LogWriter::logAndPrintError("invalid temporal_choice @U2NITS::Time::EvolveHost_2.\n");
            exit(-1);
        }
    }




    U2NITS::Space::Restrict::modifyElementFieldU2d(s->element_host, s->elementField_host);

}

void CIteration::iteration_device_20240517(myfloat& t, myfloat T, GPU::GPUSolver2* s) {
    /*
    0517GPU�����޸���
    �벻Ҫgit����Ϊ��ɫ��ʾ�޸Ĳ��֣����Կ�����Щ���޸��ˡ��ÿ��ձ��漴��
    ����ı��ˣ���ɾ��update_device_ruvp_Uold��get_global_dt_device

    �ƻ���ɾ������Ҫ�Ŀ���
    update_host_ruvp_Uold Ӧ��Ϊupdate_device_ruvp_Uold��elementField_device.ruvp�ڼ���ʱ�䲽��ʱ���õ�
    get_global_dt_host Ӧ��Ϊget_global_dt_device�����Ҫ����Ļ����͸�������dt��host(�ò��ִ������ת�Ƶ������ط�)
    ���� GPU::ElementFieldSoA::cuda_memcpy ��Ҫɾ����ת�Ƶ������ط�

    �޸���ɺ���֤��Ȼ���޸ı�������
    */

    constexpr bool debug_host_on = false;// �Ƿ���host���е��ԣ��Ա�reduce���
    if (debug_host_on) {
        //update_host_ruvp_Uold(s->elementField_host.U, s->element_vruvp, s->element_U_old, s->elementField_host.num, GlobalPara::constant::gamma);
        update_host_ruvp_Uold(s->elementField_host.U, s->elementField_host.ruvp, s->elementField_host.Uold, s->elementField_host.num, GlobalPara::constant::gamma);
    }

    //���޸�;// s->element_U_oldӦ��Ϊs->elementField_device.Uold ������elementField_device�����Uold��Ա
    update_device_ruvp_Uold(s->elementField_device, GlobalPara::constant::gamma);

    /*
    ʱ���ƽ�
    */
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
        myfloat dt = 0.0;
        if (debug_host_on) {
            // ����ģʽ����ʱ������Լ����dt��t��û�б��޸�
            dt = U2NITS::Time::get_global_dt_host(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Pr, GlobalPara::time::CFL,
                GlobalPara::constant::R, s->elementFieldVariable_dt_host, s->element_host, s->edge_host, s->elementField_host.ruvp, GlobalPara::physicsModel::equation);
        }
        else {
            myfloat dt2 = 0.0;
            // ��Լ
            GPU::Time::get_global_dt_device(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::time::CFL,
                GlobalPara::constant::R, s->elementFieldVariable_dt_device, s->element_device, s->edge_device, s->elementField_device, GlobalPara::physicsModel::equation);
            // ��Լ���
            cudaMemcpy(&dt2, s->elementFieldVariable_dt_device.dev_output, sizeof(myfloat), cudaMemcpyDeviceToHost);
            dt = dt2;
        }
        
        
        // �޸ĵ�ǰ����ʱ��t
        CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
        pCDS->set_t(t);
        pCDS->set_dt(0.0);

        t += dt;

        switch (temporal_choice) {
        case _EVO_explicit:
            GPU::Time::evolveSingleStep_device(dt, s->element_device, s->node_device, s->edge_device, s->elementField_device);
            break;
        case _EVO_rk3:
            func_unimplemented("rk3");
            break;
        default:
            func_unimplemented("temporal_choice");
            break;
        }

        if (debug_host_on) {
            // device to host���ú������ڲ��Թ�Լ������ȷ��ʱ�������������������ÿ��������
            GPU::ElementFieldSoA::cuda_memcpy(&s->elementField_host, &s->elementField_device, cudaMemcpyDeviceToHost);
        }
    }


}

void CIteration::iteration_device_old(myfloat& t, myfloat T, GPU::GPUSolver2* s) {

    update_host_ruvp_Uold(s->elementField_host.U, s->elementField_host.ruvp, s->elementField_host.Uold, s->elementField_host.num, GlobalPara::constant::gamma);

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
        myfloat dt = U2NITS::Time::get_global_dt_host(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Pr, GlobalPara::time::CFL,
            GlobalPara::constant::R, s->elementFieldVariable_dt_host, s->element_host, s->edge_host, s->elementField_host.ruvp, GlobalPara::physicsModel::equation);

        // �޸ĵ�ǰ����ʱ��t
        CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
        pCDS->set_t(t);
        pCDS->set_dt(0.0);
        t += dt;

        GPU::ElementFieldSoA::cuda_memcpy(&s->elementField_device, &s->elementField_host, cudaMemcpyHostToDevice);

        switch (temporal_choice) {
        case _EVO_explicit:
            GPU::Time::evolveSingleStep_device(dt, s->element_device, s->node_device, s->edge_device, s->elementField_device);
            break;
        case _EVO_rk3:
            func_unimplemented("rk3");
            break;
        default:
            func_unimplemented("temporal_choice");
            break;
        }

        GPU::ElementFieldSoA::cuda_memcpy(&s->elementField_host, &s->elementField_device, cudaMemcpyDeviceToHost);
    }

}

void CIteration::update_host_ruvp_Uold(const myfloat* const U[4], myfloat* ruvp[4], myfloat* U_old[4], myint num, myfloat gamma) {
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

void CIteration::update_device_ruvp_Uold(GPU::ElementFieldSoA& elementField_device, myfloat gamma) {
    GPU::Math::update_ruvp_Uold_device_2(elementField_device, gamma);
}

