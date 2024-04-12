#ifndef GPU_DATATYPE_H
#define GPU_DATATYPE_H

#ifdef _WIN32
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include <omp.h>
#endif

// define data type
// typedef作用域仅限文件，而宏(#define REAL double)作用于全局。
// 可以多次使用typedef double REAL，只要基类型都是double
// 但 typedef double REAL和typedef float REAL一起用会报错重定义
// using与typedef差不多，但赋值的方式更符合直觉；此外using可以给模板起别名

// 实数类型
using REAL = double;
using real = double;
using integer = int;

#endif // !GPU_DATATYPE_H
