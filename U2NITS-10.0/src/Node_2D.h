#ifndef NODE_2D_H
#define NODE_2D_H
#include "head.h"

class Element_2D;
class FVM_2D;

class Node_2D {
public:
	int ID = -1;
	double x = 0;
	double y = 0;
	std::vector<Element_2D*> neighborElements;//�ھӵ�Ԫ��ͨ��Ϊ3~7��

public:
	//// ����ڵ㺯��ֵ��ȡ�ھӵ�Ԫ�ķֲ�������ƽ��ֵ
	//std::vector<double> calNodeValue() const;
	// ����ڵ㺯��ֵ��ȡ�ھӵ�Ԫ�ķֲ�������ƽ��ֵ
	void calNodeValue(double res[4])const;
};

#endif