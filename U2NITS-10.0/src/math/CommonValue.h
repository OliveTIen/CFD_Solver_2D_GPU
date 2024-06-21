#ifndef MATH_COMMON_VALUE_H
#define MATH_COMMON_VALUE_H
#include "../gpu/dataType/DefineType.h"

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
		constexpr myfloat PI = 3.1415926535897;        // Բ����
		constexpr myfloat EPSILON_FLOAT_UN = 1e-44;    // float��С�ǹ�񻯸����� unnormalized
		constexpr myfloat EPSILON_FLOAT_NORMAL = 1e-37;// float��С��񻯸����� normalized
		constexpr myfloat BIG_FLOAT_NORMAL = 1e38;

		constexpr myfloat EPSILON = 1e-10;
		constexpr myfloat BIG = 1e10;

		namespace Physics {
			constexpr myfloat RHO_MAX = 300;// �ܶ����ֵ
			constexpr myfloat RHO_MIN = EPSILON;
			constexpr myfloat P_ATMOSPHERE = 101325;// ����ѹ
			constexpr myfloat P_MAX = 400 * P_ATMOSPHERE;
			constexpr myfloat P_MIN = EPSILON;

		}
	}
}

#endif 
