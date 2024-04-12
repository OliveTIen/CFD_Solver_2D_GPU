#ifndef GPU_DATATYPE_H
#define GPU_DATATYPE_H

#ifdef _WIN32
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include <omp.h>
#endif

// define data type
// typedef����������ļ�������(#define REAL double)������ȫ�֡�
// ���Զ��ʹ��typedef double REAL��ֻҪ�����Ͷ���double
// �� typedef double REAL��typedef float REALһ���ûᱨ���ض���
// using��typedef��࣬����ֵ�ķ�ʽ������ֱ��������using���Ը�ģ�������

// ʵ������
using REAL = double;
using real = double;
using integer = int;

#endif // !GPU_DATATYPE_H
