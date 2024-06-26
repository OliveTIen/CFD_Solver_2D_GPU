#include "Element_2D.h"
#include "FVM_2D.h"
#include "Edge_2D.h"
#include "output/LogWriter.h"


myfloat Element_2D::calArea(FVM_2D* f) {
    // 叉乘法计算三角形面积。取绝对值，保证结果非负
    myfloat xn[3]{}, yn[3]{};
    for (int i = 0; i < 3; i++) {
        xn[i] = f->getNodeByID(nodes[i])->x;
        yn[i] = f->getNodeByID(nodes[i])->y;
    }
    return 0.5 * abs(xn[0] * (yn[1] - yn[2]) + xn[1] * (yn[2] - yn[0]) + xn[2] * (yn[0] - yn[1]));

}


std::vector<Element_2D*> Element_2D::findNeighbor() {
    //使用前需保证pEdges是最新的，且pEdges的pElement_L/R是最新的
    std::vector<Element_2D*> pNeighborElements(3);
    for (int iedge = 0; iedge < 3; iedge++) {
        if (pEdges[iedge] == nullptr) {
            std::cout << "Error: Element_2D::pEdges[" << iedge << "] is uninitialized(Element ID: " << this->ID << ")" << std::endl;
            pNeighborElements[iedge] = nullptr;
            //continue;//结束本次循环
        }
        else {
            if (pEdges[iedge]->pElement_L == this) {
                if (pEdges[iedge]->pElement_R == nullptr) {
                    pNeighborElements[iedge] = nullptr;
                }
                else {
                    pNeighborElements[iedge] = pEdges[iedge]->pElement_R;
                }
            }
            else if(pEdges[iedge]->pElement_R == this) {
                if (pEdges[iedge]->pElement_L == nullptr) {
                    pNeighborElements[iedge] = nullptr;
                }
                else {
                    pNeighborElements[iedge] = pEdges[iedge]->pElement_L;
                }
            }
        }
    }
    return pNeighborElements;
}

std::vector<Element_2D*> Element_2D::findNeighbor_withoutNullptr() {
    std::vector<Element_2D*>neighbors_origin = this->findNeighbor();
    std::vector<Element_2D*>neighbors;
    //清除nullptr
    for (int i = 0; i < 3; i++) {
        if (neighbors_origin[i] != nullptr) {
            //neighbors_origin[i]->calxy(f);//初始化坐标 calxy!!!
            neighbors.push_back(neighbors_origin[i]);
        }
    }
    return neighbors;
}
