#ifndef MATH_COMMON_H
#define MATH_COMMON_H
#include "../gpu/dataType/Define.h"

// Common��Math��MathGPU����
namespace U2NITS {
	namespace Math {
		/*
		float(32λ����)
		�μ� csappP72 �� https://zhuanlan.zhihu.com/p/656543002
		��С�ǹ���� = 1.4e-45
		���ǹ������С����С����� �� 1.2e-38
		��С����� = 1.2e-38
		������� = 3.4e38
		*/
		constexpr REAL EPSILON_FLOAT_UN = 1e-44;// float��С�ǹ�񻯸����� unnormalized
		constexpr REAL EPSILON_FLOAT_NORMAL = 1e-37;// float��С��񻯸����� normalized
		constexpr REAL BIG_FLOAT_NORMAL = 1e38;

		constexpr REAL EPSILON = 1e-10;
		constexpr REAL BIG = 1e10;

		namespace Physics {
			constexpr real RHO_MAX = 300;// �ܶ����ֵ
			constexpr real RHO_MIN = EPSILON;
			constexpr real P_ATMOSPHERE = 101325;// ����ѹ
			constexpr real P_MAX = 400 * P_ATMOSPHERE;
			constexpr real P_MIN = EPSILON;

		}
	}
}

#endif 
