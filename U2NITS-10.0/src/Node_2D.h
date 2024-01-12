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
	std::vector<Element_T3*> neighborElements;//�ھӵ�Ԫ��ͨ��Ϊ3~7��

public:
	//����ڵ㺯��ֵ��ȡ�ھӵ�Ԫ�ķֲ�������ƽ��ֵ
	std::vector<double> calNodeValue(FVM_2D* f);
};

#endif