#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H
/**
* 线性重构
* 参照Laminar-WangQ::reconstruction_initialize.f90
*/
class FVM_2D;
class Element_2D;

class Reconstructor {
public:
	// [未使用]
	static void compute_reconstruction_matrix();
	// [过渡函数]
	static void Element_T3_updateSlope_Barth(FVM_2D* f, Element_2D* pE);

private:
	// [未完成]高斯消元法解线性方程组
	static void solveLinearEquations(double** A, double* X, double* b, int A_size_1, int A_size_2);
	static void swapRows(double* row1, double* row2, int size);
};


#endif // !RECONSTRUCTOR_H