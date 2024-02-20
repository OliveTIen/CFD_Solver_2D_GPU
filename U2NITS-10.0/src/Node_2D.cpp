#include "FVM_2D.h"

//std::vector<double> Node_2D::calNodeValue() const {
//	// 常成员函数不能修改成员变量，也不能调用非常成员函数
//	double U_node[4]{};
//	std::vector<double> uu(4, 0.0);//初始化为4个double，且赋值为0.0
//	for (int ie = 0; ie < neighborElements.size(); ie++) {
//		// 计算节点函数值。取邻居单元的分布函数的平均值
//		neighborElements[ie]->get_U(this->x, this->y, U_node);
//		for (int i = 0; i < 4; i++) {
//			uu[i] += U_node[i];
//		}
//	}
//	for (int i = 0; i < 4; i++) {
//		uu[i] /= double(neighborElements.size());
//	}
//	return uu;
//}

void Node_2D::calNodeValue(double res[4]) const {
	// 这里接受1个参数，即数组double res[4]，等效写法如下(接受两个参数)
	// double res[], int n
	// double* res, int n
	// https://blog.csdn.net/jackystudio/article/details/12137477


	double sum[4]{};// {}表示赋值为0.0 https://blog.csdn.net/ls1300005/article/details/131730596
	int nNeighbor = (int) neighborElements.size();
	for (int ie = 0; ie < nNeighbor; ie++) {
		// 获取邻居U
		double neighborElementU[4]{};
		neighborElements[ie]->get_U(this->x, this->y, neighborElementU);
		// 求和
		for (int i = 0; i < 4; i++) {
			sum[i] += neighborElementU[i];
		}
	}
	// 平均
	for (int i = 0; i < 4; i++) {
		res[i] = sum[i] / double(nNeighbor);
	}
}
