#include "../../include/drivers/CDriverBase.hpp"

CDriverBase::CDriverBase(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator) {
	rank = 0;
	size = 0;
	std::cout << "rank and size set to 0.\n";
}

CDriverBase::~CDriverBase() = default;

void CDriverBase::InitializeContainers() {

	//config_container = new CConfig * [nZone]();
	//output_container = new COutput * [nZone]();
	//geometry_container = new CGeometry * **[nZone]();
	//solver_container = new CSolver * ***[nZone]();
	//numerics_container = new CNumerics * ****[nZone]();
	//surface_movement = new CSurfaceMovement * [nZone]();
	//grid_movement = new CVolumetricMovement * *[nZone]();

	//nInst = new unsigned short[nZone];

	//for (iZone = 0; iZone < nZone; iZone++) {
	//	nInst[iZone] = 1;
	//}
}
