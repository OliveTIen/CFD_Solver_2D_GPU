#ifndef GPU_DATATYPE_H
#define GPU_DATATYPE_H

#ifdef _WIN32
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include <omp.h>
#endif

// define data type
// 相比宏定义，typedef作用域仅限定义的文件中，而宏定义作用域为全局。
typedef double REAL;// 定义实数类型为double 等价于 #define REAL double
//constexpr bool is_debugging = true; // 定义调试开关，默认为true，在编译时可通过-D_DEBUG关闭。

#endif // !GPU_DATATYPE_H
