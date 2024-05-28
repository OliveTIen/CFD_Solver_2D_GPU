#include "ResidualGPU.h"
#include "../LogWriter.h"

void GPU::Output::get_residual_device(const GPU::ElementFieldSoA& elementField_device, myfloat* residual_device, U2NITS::Output::NormType NORM_TYPE) {
	/*
	05-17 
	Ŀǰresidual��Ȼ��host�ϼ��㡣ÿ��Ҫ����в�ʱ����device���Ƶ�host��Ȼ����㡣�Ժ���Կ������ƴ˺�����ʵ��GPU����в�

	�ƻ�
	���ù�Լ����������NormType���ò�ͬ������
	����������⿪�ٿռ䣬���Ծ���flux�Ϲ�Լ����Ҫ��������n/2�Ŀռ����ڴ洢output
	��ʱӦע��֮���flux�Ͳ������ˣ���Ϊ��Ϣ����ʧ��

	norm_1 ����ֵ֮�͡������Ԫ�ؾ���ֵ����add���ӹ�Լ
	norm_2 ��Ԫ��ƽ������add���ӹ�Լ��Ȼ��
	infinity

	*/
	LogWriter::logAndPrintError("unimplemented function, GPU::Output::get_residual_device, " __FILE__);
	exit(-1);
}
