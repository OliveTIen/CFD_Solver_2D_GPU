#ifndef GPU_DATATYPE_H
#define GPU_DATATYPE_H

#ifdef _WIN32
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#endif

// define data type
// ��Ⱥ궨�壬typedef��������޶�����ļ��У����궨��������Ϊȫ�֡�
typedef double REAL;// ����ʵ������Ϊdouble �ȼ��� #define REAL double


#endif // !GPU_DATATYPE_H
