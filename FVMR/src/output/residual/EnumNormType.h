#ifndef RESIDUAL_NORM_TYPE_H
#define RESIDUAL_NORM_TYPE_H
/*
2024-05-17
为了复用NormType枚举类型，将EnumNormType专门列入一个文件
*/
namespace U2NITS {
	namespace Output {
		// 范数类型
		enum NormType {
			normType_infinity,
			normType_1,
			normType_2
		};
	}
}

#endif // !RESIDUAL_NORM_TYPE_H
