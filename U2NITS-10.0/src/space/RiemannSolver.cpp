#include "RiemannSolver.h"
#include "../global/GlobalPara.h"
#include "../Math.h"


int RiemannSolver::solve(const double* UL, const double* UR,
    const double nx, const double ny, const double length, double* flux,
    const int conservation_scheme) {
    int returnCode = -1;
    if (conservation_scheme == 
        _SOL_LocalLaxFriedrichs) {
        returnCode = LocalLaxFriedrichs(UL, UR, nx, ny, length, flux);
    }
    else {
        returnCode = ReturnStatus::invalid_solver_type;
    }
    return returnCode;
}

int RiemannSolver::LocalLaxFriedrichs(const double* UL, const double* UR,
	const double nx, const double ny, const double length, double* flux) {
    //���ܣ���ճͨ�����������������ULUR�Ȳ�������������任��������ճ��ֵͨ��flux
    //�����flux
    const double& gamma = GlobalPara::constant::gamma;
    double ruvpL[4], ruvpR[4];
    Math_2D::U_2_ruvp(UL, ruvpL, gamma);
    Math_2D::U_2_ruvp(UR, ruvpR, gamma);
    double rL = ruvpL[0];
    double uL = ruvpL[1];
    double vL = ruvpL[2];
    double pL = ruvpL[3];
    double rR = ruvpR[0];
    double uR = ruvpR[1];
    double vR = ruvpR[2];
    double pR = ruvpR[3];

    //����任 ת��Ϊ׼һά
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

    //��ճ��ֵͨ�� Roe-Pike��������� ��Ǭ����P31
    // �غ��ʽ��Local Lax-Friedrichs��ʽ
    double lambdaMax = Math::max_(abs(unL) + sqrt(gamma * pL / rL), abs(unR) + sqrt(gamma * pR / rR));
    for (int i = 0; i < 4; i++) {
        flux[i] = 0.5 * (FnL[i] + FnR[i] - lambdaMax * (UR[i] - UL[i]));
        flux[i] *= length;
    }

    return ReturnStatus::normal;
}
