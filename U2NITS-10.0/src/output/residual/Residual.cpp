#include "Residual.h"
#include "../../math/Math.h"

void U2NITS::Output::get_residual_host(const GPU::ElementFieldSoA& elementField_host, myfloat* residual_host, NormType NORM_TYPE) {
	/*
	20240411
	ֱ�ӽ��Ҷ�����Ϊ�в�
	20240517
	���Ҫ�����GPU���漰����Լ���⡣����GPU��Լ�ã�Ȼ�����ս��������host
	*/
	const myint numElements = elementField_host.num;

	for (int i = 0; i < 4; i++) {
		myfloat norm = 0;
		switch (NORM_TYPE) {
		case normType_1:// ����ֵ֮��
			for (myint j = 0; j < numElements; j++) {
				myfloat flux = U2NITS::Math::abs(elementField_host.Flux[i][j]);
				norm += flux;
			}

			break;
		case normType_2:// ƽ��֮���ٿ���
			for (myint j = 0; j < numElements; j++) {
				myfloat flux = U2NITS::Math::abs(elementField_host.Flux[i][j]);
				norm += flux * flux;
			}
			norm = sqrt(norm);

			break;
		case normType_infinity:// ���ֵ
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
