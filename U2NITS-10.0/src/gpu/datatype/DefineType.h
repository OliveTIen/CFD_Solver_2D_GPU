#ifndef DEFINE_TYPE_H
#define DEFINE_TYPE_H

// define data type
// typedef����������ļ�������(#define REAL double)������ȫ�֡�
// ���Զ��ʹ��typedef double REAL��ֻҪ�����Ͷ���double
// �� typedef double REAL��typedef float REALһ���ûᱨ���ض���
// using��typedef��࣬����ֵ�ķ�ʽ������ֱ��������using���Ը�ģ�������

// ʵ������
using REAL = float;
using myfloat = REAL;
using integer = int;
using myint = integer;

// ����ָ�룬float���͵�˫Ŀ�����
typedef void(*func_bin_myfloat) (myfloat&, const myfloat&);

#endif