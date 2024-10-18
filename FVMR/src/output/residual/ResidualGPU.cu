#include "ResidualGPU.h"
#include "../LogWriter.h"

void GPU::Output::get_residual_device(const GPU::ElementFieldSoA& elementField_device, myfloat* residual_device, U2NITS::Output::NormType NORM_TYPE) {
	/*
	05-17 
	目前residual仍然在host上计算。每次要计算残差时，从device复制到host，然后计算。以后可以考虑完善此函数，实现GPU计算残差

	计划
	调用规约函数，根据NormType设置不同的算子
	如果不想另外开辟空间，可以就在flux上规约，但要另外申请n/2的空间用于存储output
	此时应注意之后的flux就不能用了，因为信息都损失了

	norm_1 绝对值之和。先求各元素绝对值，用add算子规约
	norm_2 各元素平方，用add算子规约，然后
	infinity

	*/
	LogWriter::logAndPrintError("unimplemented function, GPU::Output::get_residual_device, " __FILE__);
	exit(-1);
}
