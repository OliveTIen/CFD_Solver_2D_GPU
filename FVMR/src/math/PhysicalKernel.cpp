#include "PhysicalKernel.h"

void U2NITS::Math::U_to_ruvp_in_batch(const myfloat* U[4], myfloat* ruvp[4], int length, myfloat gamma) {
    for (int j = 0; j < length; j++) {
        /*
        ruvp[0] = U[0];
        ruvp[1] = U[1] / U[0];
        ruvp[2] = U[2] / U[0];
        ruvp[3] = ruvp[0] * (gamma - 1) * (U[3] / U[0] - (ruvp[1] * ruvp[1] + ruvp[2] * ruvp[2]) * 0.5);
        */
        ruvp[0][j] = U[0][j];
        ruvp[1][j] = U[1][j] / U[0][j];
        ruvp[2][j] = U[2][j] / U[0][j];
        ruvp[3][j] = ruvp[0][j] * (gamma - 1) * (U[3][j] / U[0][j] - (ruvp[1][j] * ruvp[1][j] + ruvp[2][j] * ruvp[2][j]) * 0.5);
    }
}
