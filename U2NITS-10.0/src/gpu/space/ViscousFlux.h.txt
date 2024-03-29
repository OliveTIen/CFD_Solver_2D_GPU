#ifndef __VISCOUSFLUX_CUH__
#define __VISCOUSFLUX_CUH__

#include "../datatype/Define.h"
#include "../dataType/ElementSoA.h"
#include "../dataType/EdgeSoA.h"
#include "../dataType/FieldSoA.h"

namespace GPU {
	namespace Space {
		void calculateViscousFluxHost(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host);
	}
}

#endif // __VISCOUSFLUX_CUH__