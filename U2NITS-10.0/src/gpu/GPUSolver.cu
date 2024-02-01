#include "GPUSolver.h"
#include "../FVM_2D.h"
#include "GPU_space.h"


void GPUSolver::initialze() {
	/*
	将FVM_2D数据转化为Element结构并初始化
	初始化的数据有：
	    单元坐标、U、单元邻居坐标、单元邻居U
	暂时未初始化：
		单元节点坐标、边

	读取文件后，将FVM_2D的数据结构转换为Element结构：
	1. 初始化selfxy, node1, node2, node3 （由于遗留问题，self是在CPU上计算的，以后会转移到GPU上）
	2. 初始化selfValue.U1~U4, neighbor1.U1~U4, neighbor2.U1~U4, neighbor3.U1~U4
	3. 梯度初始化为0
	*/
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;
	int num_element = pFVM2D->elements.size();
	this->element_host.alloc(num_element);

	#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {
		
		Element_2D& element_i = pFVM2D->elements[i];

		// 初始化单元坐标、单元U、单元梯度
		element_host.self.isNull[i] = false;
		element_host.self.x[i] = element_i.x;
		element_host.self.y[i] = element_i.y;
		element_host.self.U1[i] = element_i.U[0];
		element_host.self.U2[i] = element_i.U[1];
		element_host.self.U3[i] = element_i.U[2];
		element_host.self.U4[i] = element_i.U[3];
		element_host.self.Ux1[i] = 0;
		element_host.self.Uy1[i] = 0;
		element_host.self.Ux2[i] = 0;
		element_host.self.Uy2[i] = 0;
		element_host.self.Ux3[i] = 0;
		element_host.self.Uy3[i] = 0;
		element_host.self.Ux4[i] = 0;
		element_host.self.Uy4[i] = 0;		

		// 初始化邻居坐标、邻居U、邻居梯度
		const int num_neighbor = 3;
		std::vector<Element_2D*> neighbors_element_i = element_i.findNeighbor();// 判断邻居是否为nullptr，并赋值。此处默认为3个邻居 ! 可能有隐患
		for (int j = 0; j < num_neighbor; j++) {
			if (neighbors_element_i[j] == nullptr) {
				element_host.neighbors[j].isNull[i] = true;
			}
			else {
				element_host.neighbors[j].isNull[i] = false;
				element_host.neighbors[j].x[i] = neighbors_element_i[j]->x;
				element_host.neighbors[j].y[i] = neighbors_element_i[j]->y;
				element_host.neighbors[j].U1[i] = neighbors_element_i[j]->U[0];
				element_host.neighbors[j].U2[i] = neighbors_element_i[j]->U[1];
				element_host.neighbors[j].U3[i] = neighbors_element_i[j]->U[2];
				element_host.neighbors[j].U4[i] = neighbors_element_i[j]->U[3];

				element_host.neighbors[j].Ux1[i] = 0;
				element_host.neighbors[j].Uy1[i] = 0;
				element_host.neighbors[j].Ux2[i] = 0;
				element_host.neighbors[j].Uy2[i] = 0;
				element_host.neighbors[j].Ux3[i] = 0;
				element_host.neighbors[j].Uy3[i] = 0;
				element_host.neighbors[j].Ux4[i] = 0;
				element_host.neighbors[j].Uy4[i] = 0;
			}
		}


	}

	/*
	cudaMalloc应放在循环外，以减少开销 参加书P364
	但是计算单元和计算边是两套不同大小的数组，因此需要申请两套
	主机复制到设备

	需要初始化：
	    单元坐标
		单元节点坐标
	
	*/
	element_device.alloc(num_element);
	

}

void GPUSolver::iteration() {
	/*
    
	1. 初始化
		单元数值通量赋为0
		根据单元U初始化单元邻居U

    2. 计算单元梯度 - GPU
	    2.1
        根据坐标、U，计算梯度
        限制器修正梯度
        检查异常值
        输入：单元坐标、单元U、单元邻居坐标、单元邻居U
        输出：单元梯度

	    2.2. 单元梯度复制到CPU，然后更新邻居梯度
		
    4. 计算边界数值通量
        计算边法向量nx ny、边坐标edgex, edgey、边长edgeLength
		    需要参数：节点坐标
        计算边左右单元在边界的U值(左右极限UL、UR)
        根据边类型，修正UL、UR
        计算黎曼问题
            需要的参数：UL UR nx ny edgeLength
            输出：边界的数值通量
    5. 计算单元数值通量
        将边界数值通量加到单元数值通量上 - 规约操作，OpenMP
		(注意每次循环需要将单元数值通量初始化为0)
    6. 显式时间推进

 
	*/

	// --- 计算单元梯度 --- 
	// 输入：单元坐标、单元U、单元邻居坐标、单元邻居U
	// 输出：单元梯度、邻居梯度
	element_host.cuda_copy_to_device(&element_device);

	int block_size = 512;
	int grid_size = (element_device.num_element + block_size - 1) / block_size;
	GPU::calculateGradient(block_size, grid_size, this->element_device);

}
