#ifndef NODE_2D_H
#define NODE_2D_H
#include "head.h"

class Element_2D;
class FVM_2D;

class Node_2D {
public:
	int ID = -1;
	myfloat x = 0;
	myfloat y = 0;
	std::vector<Element_2D*> neighborElements;//邻居单元。通常为3~7个

	int GPUID = -1; //GPU编号

public:
	//// 计算节点函数值。取邻居单元的分布函数的平均值
	//std::vector<myfloat> calNodeValue() const;
	// 计算节点函数值。取邻居单元的分布函数的平均值
	void calNodeValue(myfloat res[4])const;
};

#endif