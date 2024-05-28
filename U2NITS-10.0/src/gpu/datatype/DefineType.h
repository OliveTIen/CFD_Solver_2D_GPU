#ifndef DEFINE_TYPE_H
#define DEFINE_TYPE_H

// define data type
// typedef作用域仅限文件，而宏(#define REAL double)作用于全局。
// 可以多次使用typedef double REAL，只要基类型都是double
// 但 typedef double REAL和typedef float REAL一起用会报错重定义
// using与typedef差不多，但赋值的方式更符合直觉；此外using可以给模板起别名

// 实数类型
using REAL = float;
using myfloat = REAL;
using integer = int;
using myint = integer;

// 函数指针，float类型的双目运算符
typedef void(*func_bin_myfloat) (myfloat&, const myfloat&);

#endif