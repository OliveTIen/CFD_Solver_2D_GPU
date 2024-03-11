#ifndef _ROE_GPU_H_
#define _ROE_GPU_H_
#include "../../gpu/datatype/Define.h"

namespace GPU {
	namespace Space {
		namespace Convection {
			void ConvectRoeCommon2d(
				const REAL UL[4], const REAL UR[4], const REAL faceNormal[2], const REAL faceArea, REAL faceFlux[4], REAL gamma
			);
			void RoeDissapationTerm2d(
				REAL gamma,
				REAL ruvpL[4], REAL ruvpR[4],
				const REAL faceNormal[2], REAL faceArea,
				bool bEntropyFix, REAL KEntropyFix[3], REAL p_sensor,
				REAL drRoe[4]
			);
			inline real PShockWaveSensor() { return 0.5; }

		}
	}
}

#endif