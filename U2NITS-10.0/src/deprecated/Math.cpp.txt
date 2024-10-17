//#include "head.h"
#include "Math.h"

myfloat Math::dot(std::vector<myfloat> a, std::vector<myfloat> b) {
    if (a.size() != b.size()) {
        std::cout << "Error: a.size() != b.size(), in Math::dot\n";
        return 0.0;
    }
    myfloat ret = 0.0;
    for (int i = 0; i < a.size(); i++) {
        ret += a[i] * b[i];
    }
    return ret;
}

inline myfloat Math::vector_norm_1(std::vector<myfloat> v) {
	myfloat sum = 0.0;
	for (int i = 0; i < v.size(); i++) {
		sum += abs_(v[i]);
	}
	return sum;
}

inline myfloat Math::vector_norm_2(std::vector<myfloat> v) {
	myfloat sum = 0.0;
	for (int i = 0; i < v.size(); i++) {
		sum += v[i] * v[i];
	}
	return sqrt(sum);
}

inline myfloat Math::vector_norm_infinity(std::vector<myfloat> v) {
	myfloat m_max = 0.0;
	for (int i = 0; i < v.size(); i++) {
		m_max = max_(m_max, abs_(v[i]));
	}
	return m_max;
}

void Math::vector_times_scalar(myfloat* v, const int vlength, const myfloat s) {
	for (int i = 0; i < vlength; i++) {
		v[i] *= s;
	}
}

std::string Math::timeFormat(int dt) {
	int h = dt / 3600;//专门用int的除法
	dt -= h * 3600;
	int m = dt / 60;
	dt -= m * 60;

	std::string str;
	str += std::to_string(h);
	char szBuffer[4]{};
	sprintf_s(szBuffer, _countof(szBuffer), "%02d", m);
	str += std::string(":") + szBuffer;
	sprintf_s(szBuffer, _countof(szBuffer), "%02d", dt);
	str += std::string(":") + szBuffer;
	return str;
}

myfloat Math::RK3(myfloat x, myfloat y, myfloat h, myfloat (*f)(myfloat, myfloat)) {
	//目的：解常微分方程y'=f(x,y),例如y'=y-2x/y
	//

	//找到求积节点
	myfloat k1 = (*f)(x, y);
	myfloat k2 = (*f)(x + h * 0.5, y + 0.5 * h * k1);//这里的y是当前y
	myfloat k3 = (*f)(x + h, y - 1.0 * h * k1 + 2.0 * h * k2);
	return y + h / 6.0 * (k1 + 4 * k2 + k3);
}

void Math_2D::U_2_ruvp(const myfloat* U, myfloat* ruvp, myfloat gamma) {
	//U:	rho,rho_u,rho_v,rho_E
	//ruvp: rho,u,v,p
	ruvp[0] = U[0];
	ruvp[1] = U[1] / U[0];
	ruvp[2] = U[2] / U[0];
	

	ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
}

void Math_2D::ruvp_2_U(const myfloat* ruvp, myfloat* U, myfloat gamma) {
	//myfloat rho = U[0];
	//myfloat rho_u = U[1];//rho*u
	//myfloat rho_v = U[2];
	//myfloat rho_E = U[3];

	U[0] = ruvp[0];
	U[1] = ruvp[0] * ruvp[1];
	U[2] = ruvp[0] * ruvp[2];
	U[3] = ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]);//rhoE
}

void Math_2D::U_2_F(const myfloat* U, myfloat* F, myfloat gamma) {
	//U[0]: rho 
	//U[1]: rho_u 
	//U[2]: rho_v 
	//U[3]: rho_E 
	myfloat p = 0.5 * U[0] * (gamma - 1) * (2 * U[3] / U[0] - U[1] / U[0] * U[1] / U[0] - U[2] / U[0] * U[2] / U[0]);
	F[0] = U[1];//rho*u
	F[1] = U[1] * U[1] / U[0] + p;//rho*u^2+p
	F[2] = U[1] * U[2] / U[0];//rho*u*v
	F[3] = (U[3] + p) * U[1] / U[0];//(rho*E+p)*u

	//myfloat ret2 = gamma * p / rho;
	//if (ret2 < 0) {
	//	GlobalStatic::log("Error: gamma * p / rho < 0, at Edge\n");
	//	std::cout << "Error: gamma * p / rho < 0, at Edge\n";
	//	ret2 = 0.0;
	//}
	//return sqrt(V2) + sqrt(ret2);//gamma

}

//void Math_2D::U_2_F(const Eigen::Vector4d& U, Eigen::Vector4d& F, myfloat gamma) {
//	myfloat p = 0.5 * U[0] * (gamma - 1) * (2 * U[3] / U[0] - U[1] / U[0] * U[1] / U[0] - U[2] / U[0] * U[2] / U[0]);
//	F[0] = U[1];//rho*u
//	F[1] = U[1] * U[1] / U[0] + p;//rho*u^2+p
//	F[2] = U[1] * U[2] / U[0];//rho*u*v
//	F[3] = (U[3] + p) * U[1] / U[0];//(rho*E+p)*u
//}

//void Math_2D::U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, myfloat& lambda, myfloat gamma) {
//	myfloat rho = U[0];
//	myfloat u = U[1] / U[0];
//	myfloat v = U[2] / U[0];
//	myfloat E = U[3] / U[0];
//	Math_2D::U_2_F(U, F, gamma);
//	myfloat p = Math_2D::get_p(rho, gamma, E, u, v);
//	lambda = sqrt(u * u + v * v) + sqrt(gamma * p / rho);
//}

void Math_2D::ruvp_2_F(const myfloat* ruvp, myfloat* F, myfloat gamma) {
	F[0] = ruvp[0] * ruvp[1];
	F[1] = ruvp[0] * ruvp[1] * ruvp[1] + ruvp[3];
	F[2] = ruvp[0] * ruvp[1] * ruvp[2];
	F[3] = (ruvp[3] / (gamma - 1) + 0.5 * ruvp[0] * (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) + ruvp[3]) * ruvp[1];

	

}

void Math_2D::ruvp_2_Fn_lambda_2D(const myfloat* ruvp, myfloat* Fn, myfloat& lambda, myfloat nx, myfloat ny, myfloat gamma) {
	//计算F·n
	myfloat rho = ruvp[0];
	myfloat u = ruvp[1];
	myfloat v = ruvp[2];
	myfloat p = ruvp[3];
	myfloat un = u * nx + v * ny;
	myfloat E = get_E(ruvp, gamma);
	Fn[0] = rho * un;//rho*u*nx + rho*v*ny = rho*un
	Fn[1] = rho * u * un + p * nx;//(rho*u^2+p)*nx + (rho*u*v)*ny = rho * u * un + p * nx
	Fn[2] = rho * v * un + p + ny;//rho * v * un + p + ny
	Fn[3] = (rho * E + p) * un;//(rho * E + p) * un
	lambda = abs(un) + sqrt(gamma * p / rho);
}

void Math_2D::U_2_Fn_lambda_2D(const myfloat* U, myfloat* Fn, myfloat& lambda, myfloat nx, myfloat ny, myfloat gamma) {
	myfloat ruvp[4];
	U_2_ruvp(U, ruvp, gamma);
	ruvp_2_Fn_lambda_2D(ruvp, Fn, lambda, nx, ny, gamma);
}
