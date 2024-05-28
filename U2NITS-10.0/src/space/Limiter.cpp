#include "Limiter.h"
#include "../FVM_2D.h"

void Limiter::modifySlope_Barth(Element_2D* pE) {

    //����������������Ux,Uy�����������������������㴦��Ux,Uy�����н��ԣ�������ȫ�����н���
    //����ƫ�����½�
    FVM_2D* f = FVM_2D::getInstance();
    std::vector<Element_2D*> neighbors = pE->findNeighbor();
    myfloat UNeighborRelative[3][4]{};//�ھӺ���ֵ��������ֵ�Ĳ� i-��i���ھӣ�j-jVar
    myfloat UNeighborRelativeMax = 0;
    myfloat UNeighborRelativeMin = 0;
    for (int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i] != nullptr) {
            for (int j = 0; j < 4; j++) {
                UNeighborRelative[i][j] = neighbors[i]->U[j] - pE->U[j];
                UNeighborRelativeMax = (std::max)(UNeighborRelativeMax, UNeighborRelative[i][j]);
                UNeighborRelativeMin = (std::min)(UNeighborRelativeMin, UNeighborRelative[i][j]);
            }
        }
    }
    //���㶥��ƫ��
    myfloat UNode[3][4]{};
    myfloat UNodeRelative[3][4]{};//���㺯��ֵ��������ֵ�Ĳ�
    myfloat UNodeRelativeMax = 0;
    myfloat UNodeRelativeMin = 0;
    for (int i_node = 0; i_node < 3; i_node++) {
        Node_2D* pNode = f->pNodeTable[pE->nodes[i_node]];
        pE->get_U(pNode->x, pNode->y, UNode[i_node]);
        for (int j = 0; j < 4; j++) {
            UNodeRelative[i_node][j] = UNode[i_node][j] - pE->U[j];
            UNodeRelativeMax = (std::max)(UNodeRelativeMax, UNodeRelative[i_node][j]);
            UNodeRelativeMin = (std::min)(UNodeRelativeMin, UNodeRelative[i_node][j]);
        }
    }
    //����Ux, Uy
    myfloat scale = 1;
    if (UNodeRelativeMax > UNeighborRelativeMax) {
        scale = (std::max)(scale, UNodeRelativeMax / UNeighborRelativeMax);
    }
    if (UNodeRelativeMin < UNeighborRelativeMin) {
        scale = (std::max)(scale, UNodeRelativeMin / UNeighborRelativeMin);
    }

    for (int j = 0; j < 4; j++) {
        pE->Ux[j] /= scale;
        pE->Uy[j] /= scale;
    }

}

myfloat Limiter::BAP(myfloat* a, int n, int n_wbap) {
    myfloat sum1, sum2, t;
    sum1 = 0;
    sum2 = 0;
    // n: length of array a
    for (int i = 1; i < n; i++) {
        // ��ֹ����0
        if (abs(a[i]) < 1.0e-12) {
            return 0.0;
        }
        t = a[0] / (a[i] + 1.0e-10);
        myfloat t3 = t * t * t;
        sum1 += t3;
        sum2 += t3 * t;
    }
    return a[0] * (n_wbap + sum1) / (n_wbap + sum2);
}

myfloat Limiter::TVDlimiter(myfloat var1, myfloat var2, myfloat epsm, int lim_type) {
    // ����UNITs Space\UTIL_interpolate.F90
    myfloat ka = 1.0 / 3.0;
    myfloat result{};
    switch (lim_type) {
    case 0:// ���Բ�ֵ��������������������
    {
        result = 0.5 * ((1.0 - ka) * var2 + (1.0 + ka) * var1);
        break;
    }

    case 1:// ����smooth���������ǿ糬
    {
        myfloat smooth = (2. * var1 * var2 + epsm) / (var1 * var1 + var2 * var2 + epsm);
        result = ((1. + ka * smooth) * var1 + (1. - ka * smooth) * var2) * smooth / 2.;
        break;
    }

    case 2:// ����Vanleer���ǿ糬�߳�
    {
        myfloat r = 0.0;
        if (abs(var2) > epsm * epsm) r = var1 / var2;
        result = (r + abs(r)) / (1. + abs(r)) * var2;// v1v2���ʱΪ0��ͬ��ʱΪ2*v1*v2/(v1+v2)
        break;
    }

    case 3://����mod-smooth���,�����߳�
    {
        myfloat smooth = (2. * var1 * var2 + epsm * epsm) / (var1 * var1 + var2 * var2 + epsm * epsm);
        result = ((1. + ka * smooth) * var1 + (1. - ka * smooth) * var2) * smooth / 2. * 0.5 * abs(sgn(var1) + sgn(var2));
        // fortran��sign(a,b)=ȡǰ�����ľ���ֵ��ȡ�������ķ���
        break;
    }

    case 4://����minmod,�����������⻬����������Roe����FDS��ͨ����������ʹ�ã�������AUSM�ࡢSteger-Warming��Vanleer�����ʹ��
    {
        myfloat r = 0.0;
        if (abs(var2) > epsm * epsm) r = var1 / var2;
        myfloat minmod = min((0.5 * (r + abs(r))), 1.0);
        result = minmod * var2;
        break;
    }

    default:
        std::cout << "Error: Invalid TVDlimiter type.\n";
        break;
    }
    return result;
}
