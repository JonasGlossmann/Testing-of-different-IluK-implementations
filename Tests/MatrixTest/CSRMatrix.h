//
// Created by jonas on 29.01.2025.
//

#ifndef CSRMATRIX_H
#define CSRMATRIX_H

#endif //CSRMATRIX_H

#include "../../Matrix Management/CSRMatrix.h"
#include "../../Matrix Management/DenseVector.h"
#include "../../Tests/TestingTools/TestingTools.h"
#include <math.h>
#include <cassert>

#include "../../Preconditioners/Ilu_k_matlab_order.h"



//todo should also not be here
inline void ikj_correct_value_test() {
    //Loading control matrix with solution for first matrix of matrix filesystem
    CSRMatrix LU = CSRMatrix(R"(../Tests/TestFiles/testing_problem/lvl0/ilu0_standard_solution.txt)");

    auto [matrix, rhs, startvector]  = loadStandardTestingProblem();
    CSRMatrix A(matrix);
    std::vector<MatrixEntry> internal_LU = A.iluk(0);

    auto testmtx =  LU.getFullMatrix(false);

    for (int i = 0; i< testmtx.size(); i++) {

        assert(testmtx[i].value == internal_LU[i].value );
    }
    std::cout << "passed ikj correct value test"<< std::endl;
}

//todo should not be here
inline void ilu_matlab_residual_test() {
    //Loading control matrix with solution for first matrix of matrix filesystem
    CSRMatrix LU = CSRMatrix(R"(../Tests/TestFiles/testing_problem/lvl0/ilu0_standard_solution.txt)");

    auto [matrix, rhs, startvector]  = loadStandardTestingProblem();
    Ilu_k current_order(matrix,0);
    Ilu_k_matlab_order standard_order(matrix,0);

    DenseVector vector = current_order.apply(startvector);
    DenseVector controlVector = standard_order.apply(startvector);

    for (int i = 0; i< vector.size; i++) {
        assert(controlVector.data[i] == vector.data[i]);
    }
    std::cout << "passed ilu matlab residual test"<< std::endl;
}

//todo should be somewhere else
inline void ilu_matlab_LU_decomp_test() {
    //Loading control matrix with solution for first matrix of matrix filesystem
    CSRMatrix LU = CSRMatrix(R"(../Tests/TestFiles/testing_problem/lvl0/ilu0_standard_solution.txt)");

    auto [matrix, rhs, startvector]  = loadStandardTestingProblem();
    Ilu_k normal_order(matrix,0);
    Ilu_k_matlab_order modified_order(matrix,0);

    for(size_t i = 0;i< normal_order.U.size();i++) {
        assert(normal_order.U[i].value == modified_order.U[i]);
    }
    for(size_t i = 0;i< normal_order.L.size();i++) {
        assert(normal_order.L[i].value == modified_order.L[i]);
    }

    std::cout << "passed ilu matlab LU decomposition test"<< std::endl;
}

void CSRMatrixTest() {
    ikj_correct_value_test();
    ilu_matlab_residual_test();
    ilu_matlab_LU_decomp_test();
}