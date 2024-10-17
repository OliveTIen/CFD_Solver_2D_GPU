#ifndef LOCAL_LAX_FRIEDRICHS_H
#define LOCAL_LAX_FRIEDRICHS_H

#include "../../gpu/datatype/DefineType.h"

namespace U2NITS {
	namespace Space {
		void LocalLaxFriedrichs(const myfloat* UL, const myfloat* UR, const myfloat nx, const myfloat ny,
			const myfloat length, myfloat* flux, myfloat gamma);
	}
}

#endif // !LOCAL_LAX_FRIEDRICHS_H
