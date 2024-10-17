﻿#ifndef EDGE_2D_H
#define EDGE_2D_H
#include "head.h"

class FVM_2D;
class Element_2D;

//单元相邻的边，一维时为点
class Edge_2D {
public:
	int ID;
	int nodes[2];//节点ID
	int setID = -1;//边界ID(数组指标+1，从1开始)。通过setID可以得到setType
	Element_2D* pElement_L = nullptr;
	Element_2D* pElement_R = nullptr;
	myfloat length{};// 边长度。初始计算好，后续无需计算。为减小内存占用，用float
	myfloat refLength{};// 两侧单元中心距离，用于梯度计算。
	int GPUID = -1;

public:
	void getxy(FVM_2D* f, myfloat& x, myfloat& y);
	myfloat getx();
	myfloat gety();
	//myfloat getLength(FVM_2D* f);
	myfloat getLength();
	// 计算两个单元中心的距离
	myfloat getRefLength();

	// 未使用 计算法向，从Element_L到Element_R
	std::vector<myfloat> getDirectionN();
	// 未使用
	void getDirectionN(myfloat& nx, myfloat& ny);
	// 计算边朝向，从node0到node1
	std::vector<myfloat> getDirectionT();
	
};

#endif