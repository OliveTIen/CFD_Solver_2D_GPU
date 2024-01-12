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

double Edge_2D::getLength(FVM_2D* f) {
    double x0, y0, x1, y1;
    x0 = f->getNodeByID(nodes[0])->x;
    y0 = f->getNodeByID(nodes[0])->y;
    x1 = f->getNodeByID(nodes[1])->x;
    y1 = f->getNodeByID(nodes[1])->y;
    return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
}

double Edge_2D::getLength() {
    FVM_2D* f = FVM_2D::pFVM2D;
    double x0, y0, x1, y1;
    x0 = f->getNodeByID(nodes[0])->x;
    y0 = f->getNodeByID(nodes[0])->y;
    x1 = f->getNodeByID(nodes[1])->x;
    y1 = f->getNodeByID(nodes[1])->y;
    return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
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

std::vector<double> Edge_2D::getDirectionN(FVM_2D* f) {
    //�ɰ� ��׼ȷ
    //std::vector<double> dir(2);
    //double dx, dy, dl;
    //if (pElement_R != nullptr) {
    //    dx = pElement_R->x - pElement_L->x;
    //    dy = pElement_R->y - pElement_L->y;
    //}
    //else {
    //    double ex, ey;
    //    getxy(f, ex, ey);
    //    dx = ex - pElement_L->x;
    //    dy = ey - pElement_L->y;
    //}
    //dl = sqrt(dx * dx + dy * dy);
    //dir[0] = dx / dl;//cos(theta)
    //dir[1] = dy / dl;//sin(theta)
    //return dir;

    //�°� ֻҪgmesh���ɵ������ζ�����ʱ�붥�����У���ôע��ʱEdge��ȻҲ����ʱ�붥�����У���Element_L���ⷨ���ڱߵ��Ҳ�
    //tx=cos(theta_t),ty=sin(theta_t),
    //theta_n=theta_t-90��
    //nx=cos(theta_n)=sin=ty,ny=sin(theta_n)=-cos=-tx
    std::vector<double> T = getDirectionT(f);
    std::vector<double> N{T[1], -T[0]};
    return N;
}

std::vector<double> Edge_2D::getDirectionN() {
    std::vector<double> T = getDirectionT(FVM_2D::pFVM2D);
    std::vector<double> N{T[1], -T[0]};
    return N;
}

void Edge_2D::getDirectionN(double& nx, double& ny) {
    std::vector<double> T = getDirectionT(FVM_2D::pFVM2D);
    nx = T[1]; ny = -T[0];
}

std::vector<double> Edge_2D::getDirectionT(FVM_2D* f) {
    double dx = f->getNodeByID(nodes[1])->x - f->getNodeByID(nodes[0])->x;
    double dy = f->getNodeByID(nodes[1])->y - f->getNodeByID(nodes[0])->y;
    double dl = sqrt(dx * dx + dy * dy);
    std::vector<double> T{dx / dl, dy / dl};
    return T;
}

Eigen::Matrix4d Edge_2D::calT(FVM_2D* f,double flag) {
    Eigen::Matrix4d matT;
    std::vector<double>dir = getDirectionN(f);
    double nx = dir[0]; double ny = dir[1];
    ny *= flag;//��������sin��Ϊ�෴��
    matT <<
        1, 0, 0, 0,
        0, nx, ny, 0,
        0, -ny, nx, 0,
        0, 0, 0, 1;
    return matT;
}

Eigen::Vector4d Edge_2D::get_UL() {
    //��ȡ���е����꣬������Ԫ�ֲ���������ȡU
    double x_edge, y_edge; getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    return pElement_L->get_U(x_edge, y_edge);
}

Eigen::Vector4d Edge_2D::get_UR() {
    //��ȡ���е����꣬������Ԫ�ֲ���������ȡU
    double x_edge, y_edge; getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    return pElement_R->get_U(x_edge, y_edge);
}

void Edge_2D::get_ULUR(Eigen::Vector4d& U_L, Eigen::Vector4d& U_R) {
    //��ȡ���е����꣬������/�ҵ�Ԫ�ֲ���������ȡU
    double x_edge, y_edge; getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    U_L = pElement_L->get_U(x_edge, y_edge);
    U_R = pElement_R->get_U(x_edge, y_edge);
}

void Edge_2D::get_ULUR(double* U_L, double* U_R) {
    double x_edge, y_edge; getxy(FVM_2D::pFVM2D, x_edge, y_edge);
    pElement_L->get_U(x_edge, y_edge, U_L);
    pElement_R->get_U(x_edge, y_edge, U_R);
}
