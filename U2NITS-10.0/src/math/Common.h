#ifndef MATH_COMMON_H
#define MATH_COMMON_H
#include "../gpu/dataType/Define.h"

// Common被Math和MathGPU共用
namespace U2NITS {
	namespace Math {
		/*
		float(32位浮点)
		参见 csappP72 或 https://zhuanlan.zhihu.com/p/656543002
		最小非规格化数 = 1.4e-45
		最大非规格化数略小于最小规格化数 ≈ 1.2e-38
		最小规格化数 = 1.2e-38
		最大规格化数 = 3.4e38
		*/
		constexpr REAL EPSILON_FLOAT_UN = 1e-44;// float最小非规格化浮点数 unnormalized
		constexpr REAL EPSILON_FLOAT_NORMAL = 1e-37;// float最小规格化浮点数 normalized
		constexpr REAL BIG_FLOAT_NORMAL = 1e38;

		constexpr REAL EPSILON = 1e-10;
		constexpr REAL BIG = 1e10;

		namespace Physics {
			constexpr real RHO_MAX = 300;// 密度最大值
			constexpr real RHO_MIN = EPSILON;
			constexpr real P_ATMOSPHERE = 101325;// 大气压
			constexpr real P_MAX = 400 * P_ATMOSPHERE;
			constexpr real P_MIN = EPSILON;

		}
	}
}

#endif 
