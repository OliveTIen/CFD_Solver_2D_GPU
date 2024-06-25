#include "Evolve.h"
#include "../global/GlobalPara.h"
#include "../space/Flux.h"
#include "../output/LogWriter.h"
#include "CalculateDt.h"
#include "../space/gradient/Gradient.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"
#include "../space/viscous_flux/ViscousFluxGPU.h"



void U2NITS::Time::evolve_explicit_globaltimestep(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {

    // 计算残值。dU/dt = f(t,U)右端项
    calculateFunctionF(element_host, node_host, edge_host, elementField_host);
    // 时间推进
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < element_host.num_element; j++) {
            elementField_host.U[i][j] += dt * elementField_host.Flux[i][j];
        }
    }

}

void U2NITS::Time::evolve_rk3_globaltimestep(myfloat dt, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& yn) {
    /*
    数值求解常微分方程 dU/dt = f(t,U)，右端项即residual
    MATLAB代码如下
    function [ynp1] = rk3(xn,yn,f,h)
        % 数值求解常微分方程dy/dx = f(x,y)
        % 输入：当前x y，函数f，步长h=dx
        % 输出：下一步x
        k1=f(xn,yn);
        k2=f(xn+0.5*h,yn+0.5*h*k1);
        k3=f(xn+h, yn-h*k1+2*h*k2);
        ynp1 = yn + h/6*(k1+4*k2+k3);
    end
    事实上计算f时无需考虑xn
    以后可能会更新梯度，因此还是把FieldSoA作为一个整体复制。虽然目前梯度多复制了几份，用不上。
    如果每一小步都要更新梯度，那么就需要把前面求梯度的部分放进calculateFunctionF中
    */

    CBoundaryDoubleShockReflect* pCDS = CBoundaryDoubleShockReflect::getInstance();
    int num = element_host.num_element;
    GPU::ElementFieldSoA yn1, yn2, yn3;// k1存储于yn1.Flux，yn1存储于yn1.U
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
    calculateFunctionF(element_host, node_host, edge_host, yn1);// k1，存入yn1.Flux

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

    // 求ynp1，存入yn
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
    // 计算常微分方程的右端项f=f(t,U)。与时间无关，因此简化为f(U)

    // 根据U计算Ux、Uy，即重构
    U2NITS::Space::Gradient::Gradient(element, node, edge, ynp);
    // 根据U和Ux、Uy计算通量，存入ynp.Flux
    U2NITS::Space::Flux::calculateFluxHost(element, edge, ynp);
    // 乘以体积负倒数，得到右端项f，存入ynp.Flux
    for (int ie = 0; ie < element.num_element; ie++) {
        myfloat minus_one_on_volume = -1.0 / element.volume[ie];
        for (int j = 0; j < 4; j++) {
            ynp.Flux[j][ie] = minus_one_on_volume * ynp.Flux[j][ie];
        }
    }
}


void U2NITS::Time::evolve_explicit_localtimestep(GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host, GPU::ElementFieldVariable_dt& elementFieldVariable_dt_host) {
    // 计算残值(位于Flux)
    calculateFunctionF(element_host, node_host, edge_host, elementField_host);
    // 时间推进
    for (myint j = 0; j < element_host.num_element; j++) {
        for (int i = 0; i < 4; i++) {
            elementField_host.U[i][j] += elementFieldVariable_dt_host.alphaC[j] * elementField_host.Flux[i][j];
        }
    }

}

