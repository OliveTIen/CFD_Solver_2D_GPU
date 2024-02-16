#include "CalculateFlux.h"

void GPU::calculateFlux(ElementDataPack& element_device, EdgeSoA& edge_device) {
	// --- 计算边界数值通量 --- 
	// 输入：单元
	// 输出：边界数值通量
	// 首先初始化Flux为0，然后
	// 需要构建边界结构体，还需要初始化， 参照Solver_2D.cpp
	/*
	根据单元中心值、单元梯度、单元坐标、边坐标计算边值

	*/
}

