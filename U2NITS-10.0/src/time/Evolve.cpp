#include "Evolve.h"
#include "../global/GlobalPara.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"
#include "../math/Math.h"
#include "CalculateDt.h"
#include "../space/gradient/Gradient.h"

void U2NITS::Time::EvolveHost(REAL dt, int flag, GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
    if (flag == _EVO_explicit) {
        EvolveExplicitHost(dt, element_host, elementField_host, edge_host, edge_periodic_pair);
    }
    else {
        LogWriter::logAndPrint("Error: invalid evolve method.\n", LogWriter::Error, LogWriter::Error);
        exit(1);
    }
}

void U2NITS::Time::EvolveExplicitHost(REAL dt, GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
    // ��ֵͨ��
    U2NITS::Space::Flux::calculateFluxHost(element_host, edge_host, elementField_host);
    // ʱ���ƽ�
    for (int ie = 0; ie < element_host.num_element; ie++) {
        REAL omega = element_host.volume[ie];
        for (int j = 0; j < 4; j++) {
            elementField_host.U[j][ie] -= dt / omega * elementField_host.Flux[j][ie];
        }

    }
}

void U2NITS::Time::EvolveHost_2_addRK3(real& physicalTime, real maxPhysicalTime, GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host, GPU::NodeSoA& node_host, real* element_vruvp[4]) {
    /*
    element_vruvp ���ڼ���
    
    ���������þֲ�ʱ���ƽ���t��T��Ч
    �Ƕ�������������ʱ���ƽ�����Ҫ�ı�t��ֵ�����ȼ�����������dt��Ȼ���ƽ����ɲ���RK1��RK3
    */
    bool is_steady = GlobalPara::time::is_steady;
    bool is_explicit = GlobalPara::time::is_explicit;
    int temporal_choice = GlobalPara::time::time_advance;

    //U_to_ruvp_proMax(elementField_host.U, element_vruvp, elementField_host.num, GlobalPara::constant::gamma);

    // ����
    if (is_steady) {
        //LogWriter::logAndPrintError("unimplemented is_steady @U2NITS::Time::EvolveHost_2.\n");
        //exit(-1);
        physicalTime += 1e-3;// �ı�t��Ϊ�����tecplotʱ�ܹ������еķ�ʽ��ȡ
        evolveSteadyLocalTimeStep(element_host, node_host, edge_host, elementField_host, element_vruvp);
    }
    // �Ƕ���
    else {
        // ��������ʱ�䲽��
        REAL dt = U2NITS::Time::calculateGlobalDt(physicalTime, maxPhysicalTime, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::time::CFL, 
            GlobalPara::constant::R, element_host, edge_host, element_vruvp);
        // �޸ĵ�ǰ����ʱ��t
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

void U2NITS::Time::evolveSingleStep(real dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::FieldSoA& elementField_host) {

    int num = element_host.num_element;
    // ����dU/dt = f(t,U)�Ҷ���
    calculateFunctionF(element_host, node_host, edge_host, elementField_host);
    // ʱ�����
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < num; j++) {
            elementField_host.U[i][j] += dt * elementField_host.Flux[i][j];
        }
    }

}

void U2NITS::Time::evolveRungeKutta3(real dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::FieldSoA& yn) {
    /*
    ��ֵ��ⳣ΢�ַ��� dU/dt = f(t,U)���Ҷ��residual
    MATLAB��������
    function [ynp1] = rk3(xn,yn,f,h)
        % ��ֵ��ⳣ΢�ַ���dy/dx = f(x,y)
        % ���룺��ǰx y������f������h=dx
        % �������һ��x
        k1=f(xn,yn);
        k2=f(xn+0.5*h,yn+0.5*h*k1);
        k3=f(xn+h, yn-h*k1+2*h*k2);
        ynp1 = yn + h/6*(k1+4*k2+k3);
    end
    ��ʵ�ϼ���fʱ���迼��xn
    �Ժ���ܻ�����ݶȣ���˻��ǰ�FieldSoA��Ϊһ�����帴�ơ���ȻĿǰ�ݶȶิ���˼��ݣ��ò��ϡ�
    ���ÿһС����Ҫ�����ݶȣ���ô����Ҫ��ǰ�����ݶȵĲ��ַŽ�calculateFunctionF��
    */
    
    int num = element_host.num_element;
    GPU::FieldSoA yn1, yn2, yn3;// k1�洢��yn1.Flux��yn1�洢��yn1.U
    yn1.alloc(num);
    yn2.alloc(num);
    yn3.alloc(num);
    real* k1[4]{};
    real* k2[4]{};
    real* k3[4]{};
    for (int i = 0; i < 4; i++) {
        k1[i] = new real[num]{};
        k2[i] = new real[num]{};
        k3[i] = new real[num]{};
    }

    yn1.copyfrom(yn);
    calculateFunctionF(element_host, node_host, edge_host, yn1);// k1������yn1.Flux

    yn2.copyfrom(yn);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < num; j++) {
            yn2.U[i][j] = yn.U[i][j] + 0.5 * dt * yn1.Flux[i][j];
        }
    }
    calculateFunctionF(element_host, node_host, edge_host, yn2);// k2

    yn3.copyfrom(yn);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < num; j++) {
            yn3.U[i][j] = yn.U[i][j] - dt * yn1.Flux[i][j] + 2.0 * dt * yn2.Flux[i][j];
        }
    }
    calculateFunctionF(element_host, node_host, edge_host, yn3);// k3

    // ��ynp1������yn
    const real dt_on_6 = dt / 6.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < num; j++) {
            yn.U[i][j] = yn.U[i][j] + dt_on_6 * (yn1.Flux[i][j] + 4 * yn2.Flux[i][j] + yn3.Flux[i][j]);
        }
    }

    yn1.free();
    yn2.free();
    yn3.free();
    for (int i = 0; i < 4; i++) {
        delete[] k1[i];
        delete[] k2[i];
        delete[] k3[i];
    }
}

void U2NITS::Time::calculateFunctionF(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::FieldSoA& ynp) {
    // ���㳣΢�ַ��̵��Ҷ���f=f(t,U)����ʱ���޹أ���˼�Ϊf(U)
    
    // ����U����Ux��Uy�����ع�
    U2NITS::Space::Gradient::Gradient(element, node, edge, ynp);
    // ����U��Ux��Uy����ͨ��������ynp.Flux
    U2NITS::Space::Flux::calculateFluxHost(element, edge, ynp);
    // ����������������õ��Ҷ���f������ynp.Flux
    for (int ie = 0; ie < element.num_element; ie++) {
        real minus_one_on_volume = - 1.0 / element.volume[ie];
        for (int j = 0; j < 4; j++) {
            ynp.Flux[j][ie] = minus_one_on_volume * ynp.Flux[j][ie];
        }
    }
}


void U2NITS::Time::U_to_ruvp_proMax(const real* U[4], real* ruvp[4], int length, real gamma) {
    for (int j = 0; j < length; j++) {
        /*   
        ruvp[0] = U[0];
        ruvp[1] = U[1] / U[0];
        ruvp[2] = U[2] / U[0];
        ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
        */
        ruvp[0][j] = U[0][j];
        ruvp[1][j] = U[1][j] / U[0][j];
        ruvp[2][j] = U[2][j] / U[0][j];
        ruvp[3][j] = ruvp[0][j] * (gamma - 1) * (U[3][j] / U[0][j] - (ruvp[1][j] * ruvp[1][j] + ruvp[2][j] * ruvp[2][j]) * 0.5);
    }
}

void U2NITS::Time::evolveSteadyLocalTimeStep(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::FieldSoA& elementField, real* element_vruvp[4]) {
    /*
    �������þֲ�ʱ�䲽������������ÿ����Ԫʹ�ø��Ե�ʱ�䲽
    */
    integer num = element.num_element;
    real CFL_steady = GlobalPara::time::CFL_steady;
    calculateFunctionF(element, node, edge, elementField);
    for (int i = 0; i < 4; i++) {
        for (integer j = 0; j < num; j++) {
            real dt = 0.0;
            calculateLocalTimeStep_async_Euler(dt, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, CFL_steady,
                GlobalPara::constant::R, j, element, edge, element_vruvp);

            elementField.U[i][j] += dt * elementField.Flux[i][j];
        }
    }
}

