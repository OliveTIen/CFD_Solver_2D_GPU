/*
 * matrix.h
 * https://github.com/akalicki/matrix/blob/master/dist/matrix.h
 */

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>

class AMatrix {
    public:
        AMatrix(int, int);
        AMatrix(double**, int, int);
        AMatrix();
        ~AMatrix();
        AMatrix(const AMatrix&);
        AMatrix& operator=(const AMatrix&);

        inline double& operator()(int x, int y) { return p[x][y]; }

        AMatrix& operator+=(const AMatrix&);
        AMatrix& operator-=(const AMatrix&);
        AMatrix& operator*=(const AMatrix&);
        AMatrix& operator*=(double);
        AMatrix& operator/=(double);
        AMatrix  operator^(int);
        
        friend std::ostream& operator<<(std::ostream&, const AMatrix&);
        friend std::istream& operator>>(std::istream&, AMatrix&);

        void swapRows(int, int);
        AMatrix transpose();

        static AMatrix createIdentity(int);
        static AMatrix solve(AMatrix, AMatrix);
        static AMatrix bandSolve(AMatrix, AMatrix, int);

        // functions on vectors
        static double dotProduct(AMatrix, AMatrix);

        // functions on augmented matrices
        static AMatrix augment(AMatrix, AMatrix);
        AMatrix gaussianEliminate();
        AMatrix rowReduceFromGaussian();
        void readSolutionsFromRREF(std::ostream& os);
        AMatrix inverse();
        void print();

    private:
        int rows_, cols_;
        double **p;

        void allocSpace();
        AMatrix expHelper(const AMatrix&, int);
};

AMatrix operator+(const AMatrix&, const AMatrix&);
AMatrix operator-(const AMatrix&, const AMatrix&);
AMatrix operator*(const AMatrix&, const AMatrix&);
AMatrix operator*(const AMatrix&, double);
AMatrix operator*(double, const AMatrix&);
AMatrix operator/(const AMatrix&, double);

#endif