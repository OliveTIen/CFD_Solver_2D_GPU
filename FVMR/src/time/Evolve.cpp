#include "Evolve.h"
#include "../global/GlobalPara.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"
#include "CalculateDt.h"
#include "../space/gradient/Gradient.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"
#include "../space/viscous_flux/ViscousFluxGPU.h"



void U2NITS::Time::evolve_explicit_globaltimestep(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {

    // �����ֵ��dU/dt = f(t,U)�Ҷ���
    calculateFunctionF(element_host, node_host, edge_host, elementField_host);
    // ʱ���ƽ�
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < element_host.num_element; j++) {
            elementField_host.U[i][j] += dt * elementField_host.Flux[i][j];
        }
    }

}

void U2NITS::Time::evolve_rk3_globaltimestep(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& yn) {
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


void U2NITS::Time::evolve_explicit_localtimestep(GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_host) {
    // �����ֵ(λ��Flux)
    calculateFunctionF(element_host, node_host, edge_host, elementField_host);
    // ʱ���ƽ�
    for (myint j = 0; j < element_host.num_element; j++) {
        for (int i = 0; i < 4; i++) {
            elementField_host.U[i][j] += elementFieldVariable_dt_host.alphaC[j] * elementField_host.Flux[i][j];
        }
    }

}

