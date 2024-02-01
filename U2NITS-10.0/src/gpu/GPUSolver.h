#ifndef GPU_SOLVER_H
#define GPU_SOLVER_H

#include "DataType.h"
#include "ElementDataPack.h"

class GPUSolver {
public:
	GPU::ElementDataPack element_host;
	GPU::ElementDataPack element_device;

public:
	void initialze();

	void iteration();

    void free(){this->element_host.free();}
	//void run();

};



#endif // GPU_SOLVER_H

/*
方案1：单元仅存最少数据：
nodeCoordinates[3][3]
values[4]
neighborValues[3][3]

读取文件后，将FVM_2D的数据结构转换为Element结构：
    1. 初始化selfxy, node1, node2, node3 （由于遗留问题，self是在CPU上计算的，以后会转移到GPU上）
    2. 初始化selfValue.U1~U4, neighbor1.U1~U4, neighbor2.U1~U4, neighbor3.U1~U4
    3. 梯度初始化为0

计算梯度：（每次循环都需要计算）
    将Element传入GPU
    计算并更新单元梯度(self.Ux1~Uy4)，然后同步(barrier)，等待梯度计算完成
    导出到CPU
    在CPU中将单元梯度传递给邻居单元，更新neighbor1.Ux1~Uy4, neighbor2.Ux1~Uy4, neighbor3.Ux1~Uy4

计算边界值：
根据单元节点坐标、邻居单元梯度、邻居单元值计算边界值，同步，等待边界值计算完成
根据边界值，用黎曼求解器计算边界数值通量
根据数值通量，计算并更新单元值self.U1~U4
将单元值传递给邻居单元，更新neighbor1.U1~U4, neighbor2.U1~U4, neighbor3.U1~U4

注意：
1. Element的成员是指针，可以直接传。用到哪些就传哪些

*/