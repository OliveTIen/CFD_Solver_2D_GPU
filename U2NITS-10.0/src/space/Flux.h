#ifndef _FLUX_U2NITS_H_
#define _FLUX_U2NITS_H_

#include <map>
#include "../gpu/dataType/Datatype.h"
/*
- 曾遇到疑难问题：总是报错"不允许使用不完整的类型"，原来是ElementSoA前忘记加"GPU::"
- 20240421 非必要不暴露接口。局部函数仅定义在cpp文件中
*/

namespace U2NITS {
	namespace Space {
		namespace Flux {
			void calculateFluxHost(GPU::ElementSoA& element_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host);
		}
	}

}

#endif // !_FLUX_H_
