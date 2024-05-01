#ifndef _FLUX_U2NITS_H_
#define _FLUX_U2NITS_H_

#include <map>
#include "../gpu/dataType/Datatype.h"
/*
- �������������⣺���Ǳ���"������ʹ�ò�����������"��ԭ����ElementSoAǰ���Ǽ�"GPU::"
- 20240421 �Ǳ�Ҫ����¶�ӿڡ��ֲ�������������cpp�ļ���
*/

namespace U2NITS {
	namespace Space {
		namespace Flux {
			void calculateFluxHost(GPU::ElementSoA& element_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		}
	}

}

#endif // !_FLUX_H_
