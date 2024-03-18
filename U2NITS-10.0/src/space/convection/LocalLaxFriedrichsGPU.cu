#include "LocalLaxFriedrichsGPU.h"

__device__ void GPU::Space::Convection::LocalLaxFriedrichs2d(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny, const REAL length, REAL* flux, REAL gamma) {

    //���ܣ���ճͨ�����������������ULUR�Ȳ�������������任��������ճ��ֵͨ��flux
    //�����flux

    real ruvpL[4], ruvpR[4];
    GPU::Math::U2ruvp(UL, ruvpL, gamma);
    GPU::Math::U2ruvp(UR, ruvpR, gamma);
    real rL = ruvpL[0];
    real uL = ruvpL[1];
    real vL = ruvpL[2];
    real pL = ruvpL[3];
    real rR = ruvpR[0];
    real uR = ruvpR[1];
    real vR = ruvpR[2];
    real pR = ruvpR[3];

    //����任 ת��Ϊ׼һά
    real unL = nx * uL + ny * vL;
    real unR = nx * uR + ny * vR;

    real FnL[4]{};
    FnL[0] = rL * unL;
    FnL[1] = rL * unL * uL + nx * pL;
    FnL[2] = rL * unL * vL + ny * pL;
    real hl = pL * gamma / (rL * (gamma - 1.)) + 0.5 * (uL * uL + vL * vL);
    FnL[3] = rL * unL * hl;

    real FnR[4]{};
    FnR[0] = rR * unR;
    FnR[1] = rR * unR * uR + nx * pR;
    FnR[2] = rR * unR * vR + ny * pR;
    real hr = pR * gamma / (rR * (gamma - 1.)) + 0.5 * (uR * uR + vR * vR);
    FnR[3] = rR * unR * hr;

    //��ճ��ֵͨ�� Roe-Pike��������� ��Ǭ����P31
    // �غ��ʽ��Local Lax-Friedrichs��ʽ
    real lambdaMax = Math::max(
        abs(unL) + sqrt(gamma * pL / rL), 
        abs(unR) + sqrt(gamma * pR / rR));
    for (int i = 0; i < 4; i++) {
        flux[i] = 0.5 * (FnL[i] + FnR[i] - lambdaMax * (UR[i] - UL[i]));
        flux[i] *= length;
    }

}
