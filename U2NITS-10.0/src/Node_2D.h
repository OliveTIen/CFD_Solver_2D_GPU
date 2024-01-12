#ifndef NODE_2D_H
#define NODE_2D_H
#include "head.h"

class Element_T3;
class FVM_2D;

class Node_2D {
public:
	int ID = -1;
	double x = 0;
	double y = 0;
	std::vector<Element_T3*> neighborElements;//邻居单元。通常为3~7个

public:
	//计算节点函数值。取邻居单元的分布函数的平均值
	std::vector<double> calNodeValue(FVM_2D* f);
};

#endif