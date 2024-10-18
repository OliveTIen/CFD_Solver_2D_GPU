#include "LocalLaxFriedrichs.h"
#include "../../math/Math.h"

void U2NITS::Space::LocalLaxFriedrichs(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny, const myfloat length, myfloat* flux, myfloat gamma) {

    //功能：无粘通量黎曼求解器。根据ULUR等参数，进行坐标变换，计算无粘数值通量flux
    //输出：flux
    
    myfloat ruvpL[4], ruvpR[4];
    U2NITS::Math::U2ruvp_host(UL, ruvpL, gamma);
    U2NITS::Math::U2ruvp_host(UR, ruvpR, gamma);
    myfloat rL = ruvpL[0];
    myfloat uL = ruvpL[1];
    myfloat vL = ruvpL[2];
    myfloat pL = ruvpL[3];
    myfloat rR = ruvpR[0];
    myfloat uR = ruvpR[1];
    myfloat vR = ruvpR[2];
    myfloat pR = ruvpR[3];

    //坐标变换 转化为准一维
    myfloat unL = nx * uL + ny * vL;
    myfloat unR = nx * uR + ny * vR;

    myfloat FnL[4]{};
    FnL[0] = rL * unL;
    FnL[1] = rL * unL * uL + nx * pL;
    FnL[2] = rL * unL * vL + ny * pL;
    myfloat hl = pL * gamma / (rL * (gamma - 1.)) + 0.5 * (uL * uL + vL * vL);
    FnL[3] = rL * unL * hl;

    myfloat FnR[4]{};
    FnR[0] = rR * unR;
    FnR[1] = rR * unR * uR + nx * pR;
    FnR[2] = rR * unR * vR + ny * pR;
    myfloat hr = pR * gamma / (rR * (gamma - 1.)) + 0.5 * (uR * uR + vR * vR);
    FnR[3] = rR * unR * hr;

    //无粘数值通量 Roe-Pike黎曼求解器 王乾论文P31
    // 守恒格式：Local Lax-Friedrichs格式
    myfloat lambdaMax = (Math::max)(Math::abs(unL) + sqrt(gamma * pL / rL), Math::abs(unR) + sqrt(gamma * pR / rR));
    for (int i = 0; i < 4; i++) {
        flux[i] = 0.5 * (FnL[i] + FnR[i] - lambdaMax * (UR[i] - UL[i]));
        flux[i] *= length;
    }

}
