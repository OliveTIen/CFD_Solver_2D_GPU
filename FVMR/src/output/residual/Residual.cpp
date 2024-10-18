#include "Residual.h"
#include "../../math/Math.h"

void U2NITS::Output::get_residual_host(const GPU::ElementFieldSoA& elementField_host, myfloat* residual_host, NormType NORM_TYPE) {
	/*
	20240411
	直接将右端项作为残差
	20240517
	如果要改造成GPU，涉及到规约问题。先在GPU规约好，然后将最终结果拷贝到host
	*/
	const myint numElements = elementField_host.num;

	for (int i = 0; i < 4; i++) {
		myfloat norm = 0;
		switch (NORM_TYPE) {
		case normType_1:// 绝对值之和
			for (myint j = 0; j < numElements; j++) {
				myfloat flux = U2NITS::Math::abs(elementField_host.Flux[i][j]);
				norm += flux;
			}

			break;
		case normType_2:// 平方之和再开方
			for (myint j = 0; j < numElements; j++) {
				myfloat flux = U2NITS::Math::abs(elementField_host.Flux[i][j]);
				norm += flux * flux;
			}
			norm = sqrt(norm);

			break;
		case normType_infinity:// 最大值
			for (myint j = 0; j < numElements; j++) {
				myfloat flux = U2NITS::Math::abs(elementField_host.Flux[i][j]);
				norm = U2NITS::Math::max(norm, flux);
			}

			break;
		default:
			break;
		}
		residual_host[i] = norm;
	}
}
