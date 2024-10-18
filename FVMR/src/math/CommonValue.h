#ifndef MATH_COMMON_VALUE_H
#define MATH_COMMON_VALUE_H
#include "../gpu/dataType/DefineType.h"

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
		constexpr myfloat PI = 3.1415926535897;        // 圆周率
		constexpr myfloat EPSILON_FLOAT_UN = 1e-44;    // float最小非规格化浮点数 unnormalized
		constexpr myfloat EPSILON_FLOAT_NORMAL = 1e-37;// float最小规格化浮点数 normalized
		constexpr myfloat BIG_FLOAT_NORMAL = 1e38;

		constexpr myfloat EPSILON = 1e-10;
		constexpr myfloat BIG = 1e10;

		namespace Physics {
			constexpr myfloat RHO_MAX = 300;// 密度最大值
			constexpr myfloat RHO_MIN = EPSILON;
			constexpr myfloat P_ATMOSPHERE = 101325;// 大气压
			constexpr myfloat P_MAX = 400 * P_ATMOSPHERE;
			constexpr myfloat P_MIN = EPSILON;

		}
	}
}

#endif 
