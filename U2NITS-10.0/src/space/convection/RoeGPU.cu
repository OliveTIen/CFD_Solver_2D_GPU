#include "RoeGPU.h"

void GPU::Space::Convection::ConvectRoeCommon2d(const REAL UL[4], const REAL UR[4], const REAL faceNormal[2], const REAL faceArea, REAL faceFlux[4], REAL gamma) {
}

void GPU::Space::Convection::RoeDissapationTerm2d(REAL gamma, REAL ruvpL[4], REAL ruvpR[4], const REAL faceNormal[2], REAL faceArea, bool bEntropyFix, REAL KEntropyFix[3], REAL p_sensor, REAL drRoe[4]) {
}
