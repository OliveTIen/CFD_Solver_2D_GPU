#include "SolverDataGetter.h"
#include "GPUSolver2.h"

GPU::GPUSolver2* SolverDataGetter::getSolverInstance() {
    //auto pSolver = GPU::GPUSolver2::getInstance();
    return GPU::GPUSolver2::getInstance();
}

//GPU::OutputNodeFieldSoA* SolverDataGetter::getSolverOutputNodeFieldPointer() {
//    return &(GPU::GPUSolver2::getInstance()->outputNodeField);
//}
