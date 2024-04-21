#ifndef LOCAL_LAX_FRIEDRICHS_H
#define LOCAL_LAX_FRIEDRICHS_H

#include "../../gpu/datatype/DefineType.h"

namespace U2NITS {
	namespace Space {
		void LocalLaxFriedrichs(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny,
			const REAL length, REAL* flux, REAL gamma);
	}
}

#endif // !LOCAL_LAX_FRIEDRICHS_H
