#include "BoundaryManager.h"
#include "../global/Constexpr.h"
#include "../output/LogWriter.h"
#include "../global/GlobalPara.h"

U2NITS::BoundaryManager* U2NITS::BoundaryManager::p_instance = nullptr;

U2NITS::BoundaryManager* U2NITS::BoundaryManager::getInstance() {
	if (p_instance == nullptr) {
		p_instance = new BoundaryManager();
	}
	return p_instance;
}

int U2NITS::BoundaryManager::boundaryNameToType(std::string boundaryName) {
	// Դ�� BoundaryManager_2D::getBoundaryTypeByName
	// ����std::string�����int������GlobalPara��
	// Ӧ�����ַ���ƥ��/�ؼ��ʼ���㷨������ǰ׺����Trie������Ӧ�����Զ���ȫ��
	// �߽�����ת�߽����� �ڲ���-1
	int bType = -1;
	if (boundaryName == "wall" || boundaryName == "obstacle" || boundaryName == "airfoil" || boundaryName == "foil") {
		if (GlobalPara::physicsModel::equation == _EQ_euler)bType = _BC_wall_nonViscous;
		else bType = _BC_wall_adiabat;
	}
	else if (boundaryName == "wall_nonViscous" || boundaryName == "wall_slippery" || boundaryName == "slippery" || boundaryName == "slippery_wall") {
		bType = _BC_wall_nonViscous;
	}
	else if (boundaryName == "inlet")bType = _BC_inlet;
	else if (boundaryName == "outlet")bType = _BC_outlet;
	else if (boundaryName == "inf" || boundaryName == "infinity"
		|| boundaryName == "far" || boundaryName == "farfield")bType = _BC_inf;
	else if (boundaryName == "dsr")bType = _BC_doubleShockReflect;
	else if (boundaryName == "symmetry")bType = _BC_symmetry;
	else if (boundaryName.substr(0, 8) == "periodic") {//periodic_x
		int x = std::stoi(boundaryName.substr(boundaryName.size() - 1, 1));//ȡ���һ���ַ�
		bType = _BC_periodic_0 + x;//_BC_periodic_x = _BC_periodic_0 + x
	}
	else {
		LogWriter::logError("Error: unknown boundary type.");

		throw "Error: unknown boundary type.";
	}
	return bType;

}
