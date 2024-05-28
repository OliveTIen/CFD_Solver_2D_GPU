#include "FVM_2D.h"

//std::vector<myfloat> Node_2D::calNodeValue() const {
//	// ����Ա���������޸ĳ�Ա������Ҳ���ܵ��÷ǳ���Ա����
//	myfloat U_node[4]{};
//	std::vector<myfloat> uu(4, 0.0);//��ʼ��Ϊ4��myfloat���Ҹ�ֵΪ0.0
//	for (int ie = 0; ie < neighborElements.size(); ie++) {
//		// ����ڵ㺯��ֵ��ȡ�ھӵ�Ԫ�ķֲ�������ƽ��ֵ
//		neighborElements[ie]->get_U(this->x, this->y, U_node);
//		for (int i = 0; i < 4; i++) {
//			uu[i] += U_node[i];
//		}
//	}
//	for (int i = 0; i < 4; i++) {
//		uu[i] /= myfloat(neighborElements.size());
//	}
//	return uu;
//}

void Node_2D::calNodeValue(myfloat res[4]) const {
	// �������1��������������myfloat res[4]����Чд������(������������)
	// myfloat res[], int n
	// myfloat* res, int n
	// https://blog.csdn.net/jackystudio/article/details/12137477


	myfloat sum[4]{};// {}��ʾ��ֵΪ0.0 https://blog.csdn.net/ls1300005/article/details/131730596
	int nNeighbor = (int) neighborElements.size();
	for (int ie = 0; ie < nNeighbor; ie++) {
		// ��ȡ�ھ�U
		myfloat neighborElementU[4]{};
		neighborElements[ie]->get_U(this->x, this->y, neighborElementU);
		// ���
		for (int i = 0; i < 4; i++) {
			sum[i] += neighborElementU[i];
		}
	}
	// ƽ��
	for (int i = 0; i < 4; i++) {
		res[i] = sum[i] / myfloat(nNeighbor);
	}
}
