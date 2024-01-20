#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H
/**
* �����ع�
* ����Laminar-WangQ::reconstruction_initialize.f90
*/
class FVM_2D;
class Element_T3;

class Reconstructor {
public:
	static void compute_reconstruction_matrix();
	//[���ɺ���]
	static void Element_T3_updateSlope_Barth(FVM_2D* f, Element_T3* pE);

private:
	// [δ���]��˹��Ԫ�������Է�����
	static void solveLinearEquations(double** A, double* X, double* b, int A_size_1, int A_size_2);
	static void swapRows(double* row1, double* row2, int size);
};


#endif // !RECONSTRUCTOR_H