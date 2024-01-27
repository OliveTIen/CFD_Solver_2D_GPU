#include "Limiter.h"
#include "../FVM_2D.h"

void Limiter::modifySlope_Barth(Element_2D* pE) {

    //����������������Ux,Uy�����������������������㴦��Ux,Uy�����н��ԣ�������ȫ�����н���
    //����ƫ�����½�
    FVM_2D* f = FVM_2D::pFVM2D;
    std::vector<Element_2D*> neighbors = pE->findNeighbor();
    double UU[3][4]{};//�ھӺ���ֵ��������ֵ�Ĳ�
    double UUup = 0;
    double UUdown = 0;
    for (int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i] != nullptr) {
            for (int j = 0; j < 4; j++) {
                UU[i][j] = neighbors[i]->U[j] - pE->U[j];
                UUup = (std::max)(UUup, UU[i][j]);
                UUdown = (std::min)(UUdown, UU[i][j]);
            }
        }
    }
    //���㶥��ƫ��
    double U_node[3][4]{};
    double UU_node[3][4]{};//���㺯��ֵ��������ֵ�Ĳ�
    double UUup_node = 0;
    double UUdown_node = 0;
    for (int i_node = 0; i_node < 3; i_node++) {
        Node_2D* pNode = f->pNodeTable[pE->nodes[i_node]];
        pE->get_U(pNode->x, pNode->y, U_node[i_node]);
        for (int j = 0; j < 4; j++) {
            UU_node[i_node][j] = U_node[i_node][j] - pE->U[j];
            UUup_node = (std::max)(UUup_node, UU_node[i_node][j]);
            UUdown_node = (std::min)(UUdown_node, UU_node[i_node][j]);
        }
    }
    //����Ux, Uy
    double ratio = 1;
    if (UUup_node > UUup) {
        ratio = (std::max)(ratio, UUup_node / UUup);
    }
    if (UUdown_node < UUdown) {
        ratio = (std::max)(ratio, UUdown_node / UUdown);
    }

    for (int j = 0; j < 4; j++) {
        pE->Ux[j] /= ratio;
        pE->Uy[j] /= ratio;
    }

}

double Limiter::BAP(double* a, int n, int n_wbap) {
    double sum1, sum2, t;
    sum1 = 0;
    sum2 = 0;
    // n: length of array a
    for (int i = 1; i < n; i++) {
        // ��ֹ����0
        if (abs(a[i]) < 1.0e-12) {
            return 0.0;
        }
        t = a[0] / (a[i] + 1.0e-10);
        double t3 = t * t * t;
        sum1 += t3;
        sum2 += t3 * t;
    }
    return a[0] * (n_wbap + sum1) / (n_wbap + sum2);
}

double Limiter::TVDlimiter(double var1, double var2, double epsm, int lim_type) {
    // ����UNITs Space\UTIL_interpolate.F90
    double ka = 1.0 / 3.0;
    double result{};
    switch (lim_type) {
    case 0:// ���Բ�ֵ��������������������
    {
        result = 0.5 * ((1.0 - ka) * var2 + (1.0 + ka) * var1);
        break;
    }

    case 1:// ����smooth���������ǿ糬
    {
        double smooth = (2. * var1 * var2 + epsm) / (var1 * var1 + var2 * var2 + epsm);
        result = ((1. + ka * smooth) * var1 + (1. - ka * smooth) * var2) * smooth / 2.;
        break;
    }

    case 2:// ����Vanleer���ǿ糬�߳�
    {
        double r = 0.0;
        if (abs(var2) > epsm * epsm) r = var1 / var2;
        result = (r + abs(r)) / (1. + abs(r)) * var2;// v1v2���ʱΪ0��ͬ��ʱΪ2*v1*v2/(v1+v2)
        break;
    }

    case 3://����mod-smooth���,�����߳�
    {
        double smooth = (2. * var1 * var2 + epsm * epsm) / (var1 * var1 + var2 * var2 + epsm * epsm);
        result = ((1. + ka * smooth) * var1 + (1. - ka * smooth) * var2) * smooth / 2. * 0.5 * abs(sgn(var1) + sgn(var2));
        // fortran��sign(a,b)=ȡǰ�����ľ���ֵ��ȡ�������ķ���
        break;
    }

    case 4://����minmod,�����������⻬����������Roe����FDS��ͨ����������ʹ�ã�������AUSM�ࡢSteger-Warming��Vanleer�����ʹ��
    {
        double r = 0.0;
        if (abs(var2) > epsm * epsm) r = var1 / var2;
        double minmod = min((0.5 * (r + abs(r))), 1.0);
        result = minmod * var2;
        break;
    }

    default:
        std::cout << "Error: Invalid TVDlimiter type.\n";
        break;
    }
    return result;
}
