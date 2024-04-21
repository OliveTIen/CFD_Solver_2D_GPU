#ifndef _FLUX_U2NITS_H_
#define _FLUX_U2NITS_H_

#include <map>
#include "../gpu/dataType/Datatype.h"
/*
曾遇到疑难问题：总是报错"不允许使用不完整的类型"，原来是ElementSoA前忘记加"GPU::"

*/

namespace U2NITS {
	namespace Space {
		namespace Flux {
			using namespace GPU;
			void calculateFluxHost(ElementSoA& element_host, EdgeSoA& edge_host, FieldSoA& elementField_host);
			REAL calculateLambdaFlux(REAL edgeU[4], REAL edgeN[2], REAL gamma, REAL length);			
			void getEdgeFlux_wallNonViscous(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux);
			void getEdgeFlux_wall_adiabat(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux);
			void getEdgeFlux_farField(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* ruvp_inf, REAL* flux);
			void getEdgeFlux_farField_2(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* ruvp_inf, REAL* flux);
			void getEdgeFlux_farField_3(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, integer idx, REAL* ruvp_inf, REAL* flux);
			void modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, const REAL* ruvp_inf);
			void getEdgeFlux_inner_and_periodic(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* flux);
			void getUByXYandElementID(ElementSoA& element_host, FieldSoA& elementField_host, REAL x, REAL y, int elementID, REAL* U);
			void getEdgeFlux_farField_doubleShockReflect(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, integer idx, REAL* flux);

			void RiemannSolve(
				const REAL* UL, const REAL* UR, const REAL nx, const REAL ny,
				const REAL length, REAL* flux,
				const int scheme);

		}
	}

}

#endif // !_FLUX_H_
