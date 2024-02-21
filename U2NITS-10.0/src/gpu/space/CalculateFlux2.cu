#include "CalculateFlux2.h"

void GPU::calculateFlux2(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device) {
	// TODO: Implement this function.
}

REAL GPU::calculateLambdaFlux(REAL edgeU[4], REAL edgeN[2], REAL gamma, REAL length) {

    REAL LambdaC = 0;
    for (int ie = 0; ie < 3; ie++) {
        //eU
        auto u = edgeU[1] / edgeU[0];
        auto v = edgeU[2] / edgeU[0];
        auto E = edgeU[3] / edgeU[0];
        auto eabs = abs(u * edgeN[0] + v * edgeN[1]);//|euv¡¤en|
		//ec
		auto V2 = u * u + v * v;
		REAL& rho = edgeU[0];
		REAL p = 0.5 * rho * (gamma - 1) * (2 * E - V2);
		REAL ec = sqrt(gamma * p / rho);
        REAL& dl = length;
        LambdaC += (eabs + ec) * dl;
    }
    return LambdaC;


}
