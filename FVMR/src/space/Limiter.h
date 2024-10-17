#ifndef LIMITER_H
#define LIMITER_H
#include <iostream>
#include "../gpu/datatype/DefineType.h"
/*
[别删]
这是遗留代码。部分限制器(barth)已经转移至Gradient.cpp中
TVD还没转移
如果删了的话，以后再写TVDLimiter又要重新查一遍
*/

class Element_2D;
class Limiter {
	
public:
	static myfloat TVDlimiter(myfloat var1, myfloat var2, myfloat epsm,
		int lim_type);

private:
	static int sgn(myfloat x) {
		if (x > 0)return 1;
		else if (x == 0)return 0;
		else return -1;
	}
	static myfloat mod(myfloat a) {// abs
		return (a > 0) ? a : (-a);
	}
	static myfloat minmod(myfloat a, myfloat b) {
		if (a * b > 0)return sgn(a) * (std::min)(mod(a), mod(b));
		else return 0.0;
	}
	static myfloat vanLeer(myfloat a, myfloat b) {
		if (a == 0 && b == 0)return 0;
		return (sgn(a) + sgn(b)) * mod(a * b) / (mod(a) + mod(b));
	}

};

#endif // !LIMITER_H
