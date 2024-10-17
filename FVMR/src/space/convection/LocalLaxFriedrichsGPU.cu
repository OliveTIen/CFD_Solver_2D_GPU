#include "LocalLaxFriedrichsGPU.h"
#include "../../math/MathGPU.h"

__device__ void GPU::Space::Convection::LocalLaxFriedrichs2d(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny, const myfloat length, myfloat* flux, myfloat gamma) {

    //���ܣ���ճͨ�����������������ULUR�Ȳ�������������任��������ճ��ֵͨ��flux
    //�����flux

    myfloat ruvpL[4], ruvpR[4];
    GPU::Math::U2ruvp_device(UL, ruvpL, gamma);
    GPU::Math::U2ruvp_device(UR, ruvpR, gamma);
    myfloat rL = ruvpL[0];
    myfloat uL = ruvpL[1];
    myfloat vL = ruvpL[2];
    myfloat pL = ruvpL[3];
    myfloat rR = ruvpR[0];
    myfloat uR = ruvpR[1];
    myfloat vR = ruvpR[2];
    myfloat pR = ruvpR[3];

    //����任 ת��Ϊ׼һά
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

    //��ճ��ֵͨ�� Roe-Pike��������� ��Ǭ����P31
    // �غ��ʽ��Local Lax-Friedrichs��ʽ
    myfloat lambdaMax = Math::max(
        abs(unL) + sqrt(gamma * pL / rL), 
        abs(unR) + sqrt(gamma * pR / rR));
    for (int i = 0; i < 4; i++) {
        flux[i] = 0.5 * (FnL[i] + FnR[i] - lambdaMax * (UR[i] - UL[i]));
        flux[i] *= length;
    }

}
