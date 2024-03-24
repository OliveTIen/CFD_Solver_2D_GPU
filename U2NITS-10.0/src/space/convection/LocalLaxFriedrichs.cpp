#include "LocalLaxFriedrichs.h"
#include "../../math/Math.h"

void U2NITS::Space::LocalLaxFriedrichs(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny, const REAL length, REAL* flux, REAL gamma) {

    //功能：无粘通量黎曼求解器。根据ULUR等参数，进行坐标变换，计算无粘数值通量flux
    //输出：flux
    
    double ruvpL[4], ruvpR[4];
    U2NITS::Math::U2ruvp_host(UL, ruvpL, gamma);
    U2NITS::Math::U2ruvp_host(UR, ruvpR, gamma);
    double rL = ruvpL[0];
    double uL = ruvpL[1];
    double vL = ruvpL[2];
    double pL = ruvpL[3];
    double rR = ruvpR[0];
    double uR = ruvpR[1];
    double vR = ruvpR[2];
    double pR = ruvpR[3];

    //坐标变换 转化为准一维
    double unL = nx * uL + ny * vL;
    double unR = nx * uR + ny * vR;

    double FnL[4]{};
    FnL[0] = rL * unL;
    FnL[1] = rL * unL * uL + nx * pL;
    FnL[2] = rL * unL * vL + ny * pL;
    double hl = pL * gamma / (rL * (gamma - 1.)) + 0.5 * (uL * uL + vL * vL);
    FnL[3] = rL * unL * hl;

    double FnR[4]{};
    FnR[0] = rR * unR;
    FnR[1] = rR * unR * uR + nx * pR;
    FnR[2] = rR * unR * vR + ny * pR;
    double hr = pR * gamma / (rR * (gamma - 1.)) + 0.5 * (uR * uR + vR * vR);
    FnR[3] = rR * unR * hr;

    //无粘数值通量 Roe-Pike黎曼求解器 王乾论文P31
    // 守恒格式：Local Lax-Friedrichs格式
    double lambdaMax = (std::max)(abs(unL) + sqrt(gamma * pL / rL), abs(unR) + sqrt(gamma * pR / rR));
    for (int i = 0; i < 4; i++) {
        flux[i] = 0.5 * (FnL[i] + FnR[i] - lambdaMax * (UR[i] - UL[i]));
        flux[i] *= length;
    }

}
