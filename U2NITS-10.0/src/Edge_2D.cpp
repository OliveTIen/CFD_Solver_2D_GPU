#include "Edge_2D.h"
#include "FVM_2D.h"
#include "Math.h"

void Edge_2D::getxy(FVM_2D* f, double& x, double& y) {
    x = (f->getNodeByID(nodes[0])->x + f->getNodeByID(nodes[1])->x) / 2.0;
    y = (f->getNodeByID(nodes[0])->y + f->getNodeByID(nodes[1])->y) / 2.0;
}

double Edge_2D::getx() {
    return (FVM_2D::pFVM2D->getNodeByID(nodes[0])->x + FVM_2D::pFVM2D->getNodeByID(nodes[1])->x) / 2.0;
}

double Edge_2D::gety() {
    return (FVM_2D::pFVM2D->getNodeByID(nodes[0])->y + FVM_2D::pFVM2D->getNodeByID(nodes[1])->y) / 2.0;
}


double Edge_2D::getLength() {
    FVM_2D* f = FVM_2D::pFVM2D;
    if (f->hasInitEdgeLengths) {
        return this->length;
    }
    double x0, y0, x1, y1;
    x0 = f->getNodeByID(nodes[0])->x;
    y0 = f->getNodeByID(nodes[0])->y;
    x1 = f->getNodeByID(nodes[1])->x;
    y1 = f->getNodeByID(nodes[1])->y;
    return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
}

float Edge_2D::getRefLength() {
    // 获取两侧单元中心距离
    FVM_2D* f = FVM_2D::pFVM2D;
    if (f->hasInitEdgeLengths) {
        return this->refLength;
    }
    float xl, yl, xr, yr, dx, dy, dL;
    xl = (float)this->pElement_L->x;
    yl = (float)this->pElement_L->y;
    if (this->pElement_R == nullptr) {
        xr = (float)this->getx();
        yr = (float)this->gety();
        dx = xr - xl;
        dy = yr - yl;
        return 2.0f * sqrt(dx * dx + dy * dy);
    }
    else {
        xr = (float)this->pElement_R->x;
        yr = (float)this->pElement_R->y;
        dx = xr - xl;
        dy = yr - yl;
        return sqrt(dx * dx + dy * dy);
    }
}

void Edge_2D::U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, double& lambda) {
    double rho = U[0];
    double u = U[1] / U[0];
    double v = U[2] / U[0];
    double E = U[3] / U[0];
    double gamma = Constant::gamma;
    Math_2D::U_2_F(U, F, gamma);
    double p = Math_2D::get_p(rho, gamma, E, u, v);
    lambda = sqrt(u * u + v * v) + sqrt(gamma * p / rho);

}

std::vector<double> Edge_2D::getDirectionN() {
    std::vector<double> T = getDirectionT();
    std::vector<double> N{T[1], -T[0]};
    return N;
}

void Edge_2D::getDirectionN(double& nx, double& ny) {
    std::vector<double> T = getDirectionT();
    nx = T[1]; ny = -T[0];
}

std::vector<double> Edge_2D::getDirectionT() {
    FVM_2D* f = FVM_2D::pFVM2D;
    double dx = f->getNodeByID(nodes[1])->x - f->getNodeByID(nodes[0])->x;
    double dy = f->getNodeByID(nodes[1])->y - f->getNodeByID(nodes[0])->y;
    double dl = sqrt(dx * dx + dy * dy);
    std::vector<double> T{dx / dl, dy / dl};
    return T;
}

Eigen::Matrix4d Edge_2D::calT(FVM_2D* f,double flag) {
    Eigen::Matrix4d matT;
    std::vector<double>dir = getDirectionN();
    double nx = dir[0]; double ny = dir[1];
    ny *= flag;//逆矩阵就是sin变为相反数
    matT <<
        1, 0, 0, 0,
        0, nx, ny, 0,
        0, -ny, nx, 0,
        0, 0, 0, 1;
    return matT;
}

Eigen::Vector4d Edge_2D::get_UL() {
    //获取边中点坐标，代入左单元分布函数，获取U
    double x_edge, y_edge; getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    return pElement_L->get_U(x_edge, y_edge);
}

Eigen::Vector4d Edge_2D::get_UR() {
    //获取边中点坐标，代入左单元分布函数，获取U
    double x_edge, y_edge; getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    return pElement_R->get_U(x_edge, y_edge);
}

void Edge_2D::get_ULUR(Eigen::Vector4d& U_L, Eigen::Vector4d& U_R) {
    //获取边中点坐标，代入左/右单元分布函数，获取U
    double x_edge, y_edge; getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    U_L = pElement_L->get_U(x_edge, y_edge);
    U_R = pElement_R->get_U(x_edge, y_edge);
}

void Edge_2D::get_ULUR(double* U_L, double* U_R) {
    double x_edge, y_edge; getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    pElement_L->get_U(x_edge, y_edge, U_L);
    pElement_R->get_U(x_edge, y_edge, U_R);
}
