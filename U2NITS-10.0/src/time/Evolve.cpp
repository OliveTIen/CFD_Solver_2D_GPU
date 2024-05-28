#include "Evolve.h"
#include "../global/GlobalPara.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"
#include "../math/Math.h"
#include "CalculateDt.h"
#include "../space/gradient/Gradient.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"

void U2NITS::Time::EvolveHost_1(myfloat dt, int flag, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
    if (flag == _EVO_explicit) {
        EvolveExplicitHost(dt, element_host, elementField_host, edge_host, edge_periodic_pair);
    }
    else {
        LogWriter::logAndPrint("Error: invalid evolve method.\n", LogWriter::Error, LogWriter::Error);
        exit(1);
    }
}

void U2NITS::Time::EvolveExplicitHost(myfloat dt, GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host, std::map<int, int>& edge_periodic_pair) {
    // ��ֵͨ��
    U2NITS::Space::Flux::calculateFluxHost(element_host, edge_host, elementField_host);
    // ʱ���ƽ�
    for (int ie = 0; ie < element_host.num_element; ie++) {
        myfloat omega = element_host.volume[ie];
        for (int j = 0; j < 4; j++) {
            elementField_host.U[j][ie] -= dt / omega * elementField_host.Flux[j][ie];
        }
    }
}

void U2NITS::Time::evolve_unsteady_explicit(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {

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

void U2NITS::Time::evolve_unsteady_rk3(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& yn) {
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

    CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
    int num = element_host.num_element;
    GPU::ElementFieldSoA yn1, yn2, yn3;// k1�洢��yn1.Flux��yn1�洢��yn1.U
    yn1.alloc(num);
    yn2.alloc(num);
    yn3.alloc(num);
    myfloat* k1[4]{};
    myfloat* k2[4]{};
    myfloat* k3[4]{};
    for (int i = 0; i < 4; i++) {
        k1[i] = new myfloat[num]{};
        k2[i] = new myfloat[num]{};
        k3[i] = new myfloat[num]{};
    }

    yn1.copyfrom(yn);
    pCDS->set_dt(0.0);
    calculateFunctionF(element_host, node_host, edge_host, yn1);// k1������yn1.Flux

    yn2.copyfrom(yn);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < num; j++) {
            yn2.U[i][j] = yn.U[i][j] + 0.5 * dt * yn1.Flux[i][j];
        }
    }
    pCDS->set_dt(dt);
    calculateFunctionF(element_host, node_host, edge_host, yn2);// k2

    yn3.copyfrom(yn);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < num; j++) {
            yn3.U[i][j] = yn.U[i][j] - dt * yn1.Flux[i][j] + 2.0 * dt * yn2.Flux[i][j];
        }
    }
    pCDS->set_dt(0.5 * dt);
    calculateFunctionF(element_host, node_host, edge_host, yn3);// k3

    // ��ynp1������yn
    const myfloat dt_on_6 = dt / 6.0;
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

void U2NITS::Time::calculateFunctionF(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& ynp) {
    // ���㳣΢�ַ��̵��Ҷ���f=f(t,U)����ʱ���޹أ���˼�Ϊf(U)

    // ����U����Ux��Uy�����ع�
    U2NITS::Space::Gradient::Gradient(element, node, edge, ynp);
    // ����U��Ux��Uy����ͨ��������ynp.Flux
    U2NITS::Space::Flux::calculateFluxHost(element, edge, ynp);
    // ����������������õ��Ҷ���f������ynp.Flux
    for (int ie = 0; ie < element.num_element; ie++) {
        myfloat minus_one_on_volume = -1.0 / element.volume[ie];
        for (int j = 0; j < 4; j++) {
            ynp.Flux[j][ie] = minus_one_on_volume * ynp.Flux[j][ie];
        }
    }
}


void U2NITS::Time::U_to_ruvp_proMax(const myfloat* U[4], myfloat* ruvp[4], int length, myfloat gamma) {
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

void U2NITS::Time::evolve_steady_explicit_localTimeStep(GPU::ElementSoA& element, GPU::NodeSoA& node, GPU::EdgeSoA& edge, GPU::ElementFieldSoA& elementField, myfloat* element_vruvp[4]) {
    /*
    �������þֲ�ʱ�䲽������������ÿ����Ԫʹ�ø��Ե�ʱ�䲽
    */
    myint num = element.num_element;
    myfloat CFL_steady = GlobalPara::time::CFL_steady;
    calculateFunctionF(element, node, edge, elementField);
    for (int i = 0; i < 4; i++) {
        for (myint j = 0; j < num; j++) {
            myfloat dt = 0.0;
            calculateLocalTimeStep_async_Euler(dt, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, CFL_steady,
                GlobalPara::constant::R, j, element, edge, element_vruvp);

            elementField.U[i][j] += dt * elementField.Flux[i][j];
        }
    }
}

