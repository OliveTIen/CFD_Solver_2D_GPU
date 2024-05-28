#include "Edge_2D.h"
#include "FVM_2D.h"
#include "Math.h"

void Edge_2D::getxy(FVM_2D* f, myfloat& x, myfloat& y) {
    x = (f->getNodeByID(nodes[0])->x + f->getNodeByID(nodes[1])->x) / 2.0;
    y = (f->getNodeByID(nodes[0])->y + f->getNodeByID(nodes[1])->y) / 2.0;
}

myfloat Edge_2D::getx() {
    return (FVM_2D::getInstance()->getNodeByID(nodes[0])->x + FVM_2D::getInstance()->getNodeByID(nodes[1])->x) / 2.0;
}

myfloat Edge_2D::gety() {
    return (FVM_2D::getInstance()->getNodeByID(nodes[0])->y + FVM_2D::getInstance()->getNodeByID(nodes[1])->y) / 2.0;
}


myfloat Edge_2D::getLength() {
    FVM_2D* f = FVM_2D::getInstance();
    if (f->hasInitEdgeLengths) {
        return this->length;
    }
    myfloat x0, y0, x1, y1;
    x0 = f->getNodeByID(nodes[0])->x;
    y0 = f->getNodeByID(nodes[0])->y;
    x1 = f->getNodeByID(nodes[1])->x;
    y1 = f->getNodeByID(nodes[1])->y;
    return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
}

myfloat Edge_2D::getRefLength() {
    // 获取两侧单元中心距离
    FVM_2D* f = FVM_2D::getInstance();
    if (f->hasInitEdgeLengths) {
        return this->refLength;
    }
    float xl, yl, xr, yr, dx, dy;
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

//void Edge_2D::U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, myfloat& lambda) {
//    myfloat rho = U[0];
//    myfloat u = U[1] / U[0];
//    myfloat v = U[2] / U[0];
//    myfloat E = U[3] / U[0];
//    myfloat gamma = GlobalPara::constant::gamma;
//    Math_2D::U_2_F(U, F, gamma);
//    myfloat p = Math_2D::get_p(rho, gamma, E, u, v);
//    lambda = sqrt(u * u + v * v) + sqrt(gamma * p / rho);
//
//}

std::vector<myfloat> Edge_2D::getDirectionN() {
    std::vector<myfloat> T = getDirectionT();
    std::vector<myfloat> N{T[1], -T[0]};
    return N;
}

void Edge_2D::getDirectionN(myfloat& nx, myfloat& ny) {
    std::vector<myfloat> T = getDirectionT();
    nx = T[1]; ny = -T[0];
}

std::vector<myfloat> Edge_2D::getDirectionT() {
    FVM_2D* f = FVM_2D::getInstance();
    myfloat dx = f->getNodeByID(nodes[1])->x - f->getNodeByID(nodes[0])->x;
    myfloat dy = f->getNodeByID(nodes[1])->y - f->getNodeByID(nodes[0])->y;
    myfloat dl = sqrt(dx * dx + dy * dy);
    std::vector<myfloat> T{dx / dl, dy / dl};
    return T;
}

Eigen::Matrix4d Edge_2D::calT(FVM_2D* f,myfloat flag) {
    Eigen::Matrix4d matT;
    std::vector<myfloat>dir = getDirectionN();
    myfloat nx = dir[0]; myfloat ny = dir[1];
    ny *= flag;//逆矩阵就是sin变为相反数
    matT <<
        1, 0, 0, 0,
        0, nx, ny, 0,
        0, -ny, nx, 0,
        0, 0, 0, 1;
    return matT;
}

//Eigen::Vector4d Edge_2D::get_UL() {
//    //获取边中点坐标，代入左单元分布函数，获取U
//    myfloat x_edge, y_edge; getxy(FVM_2D::getInstance(), x_edge, y_edge);
//    return pElement_L->get_U(x_edge, y_edge);
//}

//Eigen::Vector4d Edge_2D::get_UR() {
//    //获取边中点坐标，代入左单元分布函数，获取U
//    myfloat x_edge, y_edge; getxy(FVM_2D::getInstance(), x_edge, y_edge);
//    return pElement_R->get_U(x_edge, y_edge);
//}

//void Edge_2D::get_ULUR(Eigen::Vector4d& U_L, Eigen::Vector4d& U_R) {
//    //获取边中点坐标，代入左/右单元分布函数，获取U
//    myfloat x_edge, y_edge; getxy(FVM_2D::getInstance(), x_edge, y_edge);
//    U_L = pElement_L->get_U(x_edge, y_edge);
//    U_R = pElement_R->get_U(x_edge, y_edge);
//}

void Edge_2D::get_ULUR(myfloat* U_L, myfloat* U_R) {
    myfloat x_edge, y_edge; getxy(FVM_2D::getInstance(), x_edge, y_edge);
    pElement_L->get_U(x_edge, y_edge, U_L);
    pElement_R->get_U(x_edge, y_edge, U_R);
}
