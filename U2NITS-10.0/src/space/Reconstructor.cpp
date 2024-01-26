#include "Reconstructor.h"
#include "../FVM_2D.h"
#include "../math/AMatrix.h"
#include "Limiter.h"

void Reconstructor::compute_reconstruction_matrix() {



}

void Reconstructor::solveLinearEquations(double** A, double* X, double* b, int A_size_1, int A_size_2) {
    // ��˹��Ԫ���� https://blog.csdn.net/y6123236/article/details/128476344
    // ��ʼ���������ע�⣬��matA���޸Ĳ���Ӱ��A��b��Ԫ��
    int size1 = A_size_1;
    int size2 = A_size_2 + 1;
    double** Ab = new double* [size1];
    for (int i = 0; i < size1; i++) {
        Ab[i] = new double[size2];
        for (int j = 0; j < size2 - 1; j++) {
            Ab[i][j] = A[i][j];
        }
        // �����������һ�о���b
        Ab[i][size2 - 1] = b[i];
    }
    
    // �ҳ������(��������Ԫ��Ϊ��׼)�����һ�н���
    double maxRow_num = Ab[0][0];
    int maxRow_index = 0;
    for (int i = 0; i < size1; i++) {
        if (Ab[i][0]> maxRow_num) {
            maxRow_num = Ab[i][0];
            maxRow_index = i;
        }
    }
    if (maxRow_index != 0) {
        // ��������
        swapRows(Ab[0], Ab[maxRow_index], size2);
    }



    // �ͷ�new�����ڴ�
    for (int i = 0; i < size1; ++i) {
        delete[] Ab[i];
    }
    delete[]Ab;
}

void Reconstructor::swapRows(double* row1, double* row2, int size) {
    double tmp;
    for (int j = 0; j < size; j++) {
        tmp = row1[j];
        row1[j] = row2[j];
        row2[j] = tmp;
    }
}


void Reconstructor::Element_T3_updateSlope_Barth(FVM_2D* f, Element_T3* pE) {
    //����UxUy����������
    // ÿ���غ�����б�����໥�����ģ����Բ���

    //3���ھ� ������ ��С���˷�
    //   |x0-x y0-y|        |U0i-Ui|
    //   |x1-x y1-y| |dUdx|=|U1i-Ui|   �� dX * dUdX = dU
    //   |x2-x y2-y| |dUdy| |U2i-Ui|
    //   

    //2���ھ� �����Է�����
    //   
    //   |x0-x y0-y| |Uxi|=|U0i-Ui|
    //   |x1-x y1-y| |Uyi| |U1i-Ui|
    //   

    //1���ھ� ���� ������������ֱ�����ݶ�Ϊ0
    //   
    //   |x0-x y0-y| |Uxi|=|U0i-Ui|
    //   |y-y0 x0-x| |Uyi|=|  0   |
    //   

    std::vector<Element_T3*>neighbors = pE->findNeighbor_withoutNullptr();
    const int nNeighbor = (int)neighbors.size();
    const int nVar = 4;// �غ�����������άΪ4
    double dX[3][2]{};
    double dUdX[2][nVar]{};
    double dU[3][nVar]{};
    // ��ʼ��dX��dU
    for (int j = 0; j < nNeighbor; j++) {
        dX[j][0] = neighbors[j]->x - pE->x;
        dX[j][1] = neighbors[j]->y - pE->y;
        dU[j][0] = neighbors[j]->U[0] - pE->U[0];
        dU[j][1] = neighbors[j]->U[1] - pE->U[1];
        dU[j][2] = neighbors[j]->U[2] - pE->U[2];
        dU[j][3] = neighbors[j]->U[3] - pE->U[3];
    }
    // ���dUdX
    if (nNeighbor == 3) {
        // ��С���˷� x=(A'A)^{-1}A'b ����A=dXm,b=dUm
        for (int iVar = 0; iVar < nVar; iVar++) {
            AMatrix dXm(3, 2), dUdXm(2, 1), dUm(3, 1);
            dXm(0, 0) = dX[0][0];
            dXm(0, 1) = dX[0][1];
            dXm(1, 0) = dX[1][0];
            dXm(1, 1) = dX[1][1];
            dXm(2, 0) = dX[2][0];
            dXm(2, 1) = dX[2][1];
            dUm(0, 0) = dU[0][iVar];
            dUm(1, 0) = dU[1][iVar];
            
            dUdX[0][iVar] = dUdXm(0, 0);
            dUdX[1][iVar] = dUdXm(1, 0);
            dUdXm = (dXm.transpose() * dXm).inverse() * dXm.transpose() * dUm;
            dUdX[0][iVar] = dUdXm(0, 0);
            dUdX[1][iVar] = dUdXm(1, 0);
        }

    }
    else if (nNeighbor == 2) {
        // �����Է�����
        try{
            for (int iVar = 0; iVar < nVar; iVar++) {
                AMatrix dXm(2, 2), dUdXm(2, 1), dUm(2, 1);
                dXm(0, 0) = dX[0][0];
                dXm(0, 1) = dX[0][1];
                dXm(1, 0) = dX[1][0];
                dXm(1, 1) = dX[1][1];
                dUm(0, 0) = dU[0][iVar];
                dUm(1, 0) = dU[1][iVar];
                dUdXm = AMatrix::solve(dXm, dUm);
                dUdX[0][iVar] = dUdXm(0, 0);
                dUdX[1][iVar] = dUdXm(1, 0);
            }
        }
        catch (std::domain_error e) {
            std::stringstream error_string;
            error_string << e.what() << "; ";
            error_string << "in element : ID = " << pE->ID
                << ", x=" << pE->x
                << ", y=" << pE->y
                << "\n";
            std::cout << error_string.str();
            //throw e;
        }

    }
    else if (nNeighbor == 1) {
        dX[1][0] = -dX[0][1];
        dX[1][1] = dX[0][0];
        dU[1][0] = 0;
        for (int iVar = 0; iVar < nVar; iVar++) {
            AMatrix dXm(2, 2), dUdXm(2, 1), dUm(2, 1);
            dXm(0, 0) = dX[0][0];
            dXm(0, 1) = dX[0][1];
            dXm(1, 0) = dX[1][0];
            dXm(1, 1) = dX[1][1];
            dUm(0, 0) = dU[0][iVar];
            dUm(1, 0) = dU[1][iVar];
            dUdXm = AMatrix::solve(dXm, dUm);
            dUdX[0][iVar] = dUdXm(0, 0);
            dUdX[1][iVar] = dUdXm(1, 0);
        }

    }
    else {
        std::cout << "Error: invalid neighbor number.\n";
    }
    // ��dUdX�����Ԫ
    for (int i = 0; i < nVar; i++) {
        pE->Ux[i] = dUdX[0][i];
        pE->Uy[i] = dUdX[1][i];
    }

    //Barth������
    //pE->restructor_in_updateSlope_Barth(f);
    Limiter::modifySlope_Barth(pE);

    // ����쳣ֵ
    std::stringstream error_string;
    for (int i = 0; i < nVar; i++) {
        if (isnan(pE->Ux[i])) {
            error_string << "isnan(pElement->Ux[" << i << "]), in element : ID = " << pE->ID
                << ", x=" << pE->x
                << ", y=" << pE->y
                << "\n";
        }
        if (isnan(pE->Uy[i])) {
            error_string << "isnan(pElement->Uy[" << i << "]), in element : ID = " << pE->ID
                << ", x=" << pE->x
                << ", y=" << pE->y
                << "\n";
        }
    }
    std::cout << error_string.str();

}
