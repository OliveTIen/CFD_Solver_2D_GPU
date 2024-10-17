#include "Limiter.h"
#include "../legacy/FVM_2D.h"


myfloat Limiter::TVDlimiter(myfloat var1, myfloat var2, myfloat epsm, int lim_type) {
    // 参照UNITs Space\UTIL_interpolate.F90
    myfloat ka = 1.0 / 3.0;
    myfloat result{};
    switch (lim_type) {
    case 0:// 线性插值，无限制器，低速流动
    {
        result = 0.5 * ((1.0 - ka) * var2 + (1.0 + ka) * var1);
        break;
    }

    case 1:// 二阶smooth限制器，亚跨超
    {
        myfloat smooth = (2. * var1 * var2 + epsm) / (var1 * var1 + var2 * var2 + epsm);
        result = ((1. + ka * smooth) * var1 + (1. - ka * smooth) * var2) * smooth / 2.;
        break;
    }

    case 2:// 二阶Vanleer，亚跨超高超
    {
        myfloat r = 0.0;
        if (abs(var2) > epsm * epsm) r = var1 / var2;
        result = (r + abs(r)) / (1. + abs(r)) * var2;// v1v2异号时为0，同号时为2*v1*v2/(v1+v2)
        break;
    }

    case 3://二阶mod-smooth混合,超、高超
    {
        myfloat smooth = (2. * var1 * var2 + epsm * epsm) / (var1 * var1 + var2 * var2 + epsm * epsm);
        result = ((1. + ka * smooth) * var1 + (1. - ka * smooth) * var2) * smooth / 2. * 0.5 * abs(sgn(var1) + sgn(var2));
        // fortran中sign(a,b)=取前面数的绝对值，取后面数的符号
        break;
    }

    case 4://二阶minmod,该限制器不光滑，不建议与Roe这种FDS类通量求解器配合使用，建议与AUSM类、Steger-Warming、Vanleer等配合使用
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
