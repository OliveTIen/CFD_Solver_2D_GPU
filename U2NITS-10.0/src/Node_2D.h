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
	std::vector<Element_2D*> neighborElements;//�ھӵ�Ԫ��ͨ��Ϊ3~7��

	int GPUID = -1; //GPU���

public:
	//// ����ڵ㺯��ֵ��ȡ�ھӵ�Ԫ�ķֲ�������ƽ��ֵ
	//std::vector<myfloat> calNodeValue() const;
	// ����ڵ㺯��ֵ��ȡ�ھӵ�Ԫ�ķֲ�������ƽ��ֵ
	void calNodeValue(myfloat res[4])const;
};

#endif