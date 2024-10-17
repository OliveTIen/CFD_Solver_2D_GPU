#include "Limiter.h"
#include "../legacy/FVM_2D.h"


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
