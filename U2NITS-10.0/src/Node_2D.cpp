#include "FVM_2D.h"

std::vector<double> Node_2D::calNodeValue(FVM_2D* f) {
	double U_node[4]{};
	std::vector<double> uu(4, 0.0);//初始化为4个double，且赋值为0.0
	
	for (int ie = 0; ie < neighborElements.size(); ie++) {
		neighborElements[ie]->get_U(this->x, this->y, U_node);
		for (int i = 0; i < 4; i++) {
			uu[i] += U_node[i];
		}
	}
	for (int i = 0; i < 4; i++) {
		uu[i] /= double(neighborElements.size());
	}
	return uu;
}
