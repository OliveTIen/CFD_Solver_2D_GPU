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
    学习const和指针修饰符的使用。详见README(05-17)
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
    时间推进
    定常，采用局部时间推进，t、T无效
    非定常，采用整体时间推进，需要改变t的值。首先计算允许的最大dt，然后推进。可采用RK1或RK3
    */
    using namespace U2NITS::Time;
    bool is_steady = GlobalPara::time::is_steady;
    bool is_explicit = GlobalPara::time::is_explicit;
    int temporal_choice = GlobalPara::time::time_advance;

    // 定常
    if (is_steady) {
        t += 1e-3;// 改变t是为了输出tecplot时能够以序列的方式读取，方便观察动画以进行调试
        evolve_steady_explicit_localTimeStep(s->element_host, s->node_host, s->edge_host, s->elementField_host, s->elementField_host.ruvp);
    }
    // 非定常
    else {
        // 计算整体时间步长
        myfloat dt = U2NITS::Time::get_global_dt_host(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Pr, GlobalPara::time::CFL,
            GlobalPara::constant::R, s->elementFieldVariable_dt_host, s->element_host, s->edge_host, s->elementField_host.ruvp, GlobalPara::physicsModel::equation);
        // 双马赫反射修改边界
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
    0517GPU代码修改中
    请不要git，因为蓝色表示修改部分，可以看出那些被修改了。用快照保存即可
    如果改崩了，请删除update_device_ruvp_Uold和get_global_dt_device

    计划：删除不必要的拷贝
    update_host_ruvp_Uold 应改为update_device_ruvp_Uold，elementField_device.ruvp在计算时间步长时会用到
    get_global_dt_host 应改为get_global_dt_device，如果要输出的话，就复制最终dt到host(该部分代码可以转移到其他地方)
    两处 GPU::ElementFieldSoA::cuda_memcpy 都要删除，转移到其他地方

    修改完成后，验证，然后修改本函数名
    */

    constexpr bool debug_host_on = false;// 是否开启host进行调试，对比reduce结果
    if (debug_host_on) {
        //update_host_ruvp_Uold(s->elementField_host.U, s->element_vruvp, s->element_U_old, s->elementField_host.num, GlobalPara::constant::gamma);
        update_host_ruvp_Uold(s->elementField_host.U, s->elementField_host.ruvp, s->elementField_host.Uold, s->elementField_host.num, GlobalPara::constant::gamma);
    }

    //待修改;// s->element_U_old应改为s->elementField_device.Uold 。请在elementField_device中添加Uold成员
    update_device_ruvp_Uold(s->elementField_device, GlobalPara::constant::gamma);

    /*
    时间推进
    */
    bool is_steady = GlobalPara::time::is_steady;
    bool is_explicit = GlobalPara::time::is_explicit;
    int temporal_choice = GlobalPara::time::time_advance;
    auto func_unimplemented = [](std::string str = "")->void {
        LogWriter::logAndPrintError("unimplemented " + str + " @CIteration::m_EvolveGPU_2\n");
        exit(-1);
    };

    // 定常
    if (is_steady) {
        func_unimplemented("steady");
    }
    // 非定常
    else {
        // 计算整体时间步长
        myfloat dt = 0.0;
        if (debug_host_on) {
            // 调试模式开启时，仅规约计算dt，t并没有被修改
            dt = U2NITS::Time::get_global_dt_host(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Pr, GlobalPara::time::CFL,
                GlobalPara::constant::R, s->elementFieldVariable_dt_host, s->element_host, s->edge_host, s->elementField_host.ruvp, GlobalPara::physicsModel::equation);
        }
        else {
            myfloat dt2 = 0.0;
            // 规约
            GPU::Time::get_global_dt_device(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Re, GlobalPara::constant::Pr, GlobalPara::time::CFL,
                GlobalPara::constant::R, s->elementFieldVariable_dt_device, s->element_device, s->edge_device, s->elementField_device, GlobalPara::physicsModel::equation);
            // 规约结果
            cudaMemcpy(&dt2, s->elementFieldVariable_dt_device.dev_output, sizeof(myfloat), cudaMemcpyDeviceToHost);
            dt = dt2;
        }
        
        
        // 修改当前物理时间t
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
            // device to host。该函数仅在测试规约函数正确性时开启。正常情况下无需每步都拷贝
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

    // 定常
    if (is_steady) {
        func_unimplemented("steady");
    }
    // 非定常
    else {
        // 计算整体时间步长
        myfloat dt = U2NITS::Time::get_global_dt_host(t, T, GlobalPara::constant::gamma, GlobalPara::constant::Pr, GlobalPara::time::CFL,
            GlobalPara::constant::R, s->elementFieldVariable_dt_host, s->element_host, s->edge_host, s->elementField_host.ruvp, GlobalPara::physicsModel::equation);

        // 修改当前物理时间t
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
    用U更新element_vruvp.element_vruvp用于计算时间步长Dt以及输出流场
    U:rho,rho_u,rho_v,rho_E ruvp:rho,u,v,p
    */
    for (int i = 0; i < num; i++) {
        // 更新host端的U_old
        for (int j = 0; j < 4; j++) {
            U_old[j][i] = U[j][i];
        }
        // 更新host端的ruvp  
        ruvp[0][i] = U[0][i];
        ruvp[1][i] = U[1][i] / U[0][i];
        ruvp[2][i] = U[2][i] / U[0][i];
        ruvp[3][i] = ruvp[0][i] * (gamma - 1) * (U[3][i] / U[0][i] - (ruvp[1][i] * ruvp[1][i] + ruvp[2][i] * ruvp[2][i]) * 0.5);
    }

}

void CIteration::update_device_ruvp_Uold(GPU::ElementFieldSoA& elementField_device, myfloat gamma) {
    GPU::Math::update_ruvp_Uold_device_2(elementField_device, gamma);
}

