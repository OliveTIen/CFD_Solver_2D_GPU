#ifndef LIMITER_H
#define LIMITER_H
#include <iostream>
class Element_T3;
class Limiter {
	
public:
	static void modifySlope_Barth(Element_T3* pE);
	static double BAP(double* a, int n, int n_wbap);
	static double TVDlimiter(double var1, double var2, double epsm,
		int lim_type);

private:
	static int sgn(double x) {
		if (x > 0)return 1;
		else if (x == 0)return 0;
		else return -1;
	}
	static double mod(double a) {// abs
		return (a > 0) ? a : (-a);
	}
	static double minmod(double a, double b) {
		if (a * b > 0)return sgn(a) * (std::min)(mod(a), mod(b));
		else return 0.0;
	}
	static double vanLeer(double a, double b) {
		if (a == 0 && b == 0)return 0;
		return (sgn(a) + sgn(b)) * mod(a * b) / (mod(a) + mod(b));
	}

};

#endif // !LIMITER_H
