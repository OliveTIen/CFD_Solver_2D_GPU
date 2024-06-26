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

