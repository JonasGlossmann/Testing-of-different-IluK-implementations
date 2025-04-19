//
// Created by jonas on 03.04.2025.
//

#ifndef ILUK_MODERNCPP_STANDARD_ORDER_TEST_H
#define ILUK_MODERNCPP_STANDARD_ORDER_TEST_H

#include "../../Preconditioners/Iluk_moderncpp_standard_order.h"
#include "../../Preconditioners/Ilu_k_matlab_order.h"
#include "../TestingTools/TestingTools.h"


inline void testing_grounds() {
}

/*
inline void testing_sparse_vector_subtraction() {
    std::vector<MatrixEntry> x = {{0,0,1.0}, {0,4,1.0}};
    std::vector<MatrixEntry> y = {{0,0,1.0}, {0,3,1.0}};
    std::vector<int> x_lvl = {0,0};
    std::vector<int> y_lvl = {1,1};
    double alpha = 0.5;
    int alpha_lvl = 2;

    sparse_row_subtraction_with_level(x,y,x_lvl,y_lvl,alpha,alpha_lvl);

    debug_assert(0.5, x[0].value,"expected 0.5, got: " + std::to_string(x[0].value));
    debug_assert(0,x[0].col,"expected 0, got " + std::to_string(x[0].col));
    debug_assert(-0.5, x[1].value,"expected -0.5, got " + std::to_string(x[1].value));
    debug_assert(3,x[1].col,"expected 3, got " + std::to_string(x[1].col));
    debug_assert(1,x[2].value,"expected 1, got " + std::to_string(x[2].value));
    debug_assert(4,x[2].col,"expected 4, got " + std::to_string(x[2].col));

    debug_assert(0, x_lvl[0],"expected 1, got " + std::to_string(x_lvl[0]));
    debug_assert(4, x_lvl[1],"expected 2, got " + std::to_string(x_lvl[1]));
    debug_assert(0, x_lvl[2],"expected 0, got " + std::to_string(x_lvl[2]));

    debug_assert(0,x[0].row,"expected 0, got " + std::to_string(x[0].row));
    debug_assert(0,x[1].row,"expected 0, got " + std::to_string(x[1].row));
    debug_assert(0,x[2].row,"expected 0, got " + std::to_string(x[2].row));
    std::cout << "testing_sparse_vector_subtraction"<< std::endl;
}

inline void testing_filter_entries() {
    std::vector<MatrixEntry> x = {{0,0,1.0},{0,3,2.0}, {0,4,3.0}};
    std::vector<int> x_lvl = {0,1,2};

    filter_entries(x,x_lvl,1);

    debug_assert(1, x[0].value,"expected 1, got: " + std::to_string(x[0].value));
    debug_assert(2, x[1].value,"expected 2, got " + std::to_string(x[1].value));
    debug_assert(2, x.size(),"expected 2, got " + std::to_string(x.size()));

    debug_assert(2, x_lvl.size(),"expected 2, got " + std::to_string(x_lvl.size()));
    debug_assert(0,x_lvl[0],"expected 0, got " + std::to_string(x_lvl[0]));
    debug_assert(1,x_lvl[1],"expected 1, got " + std::to_string(x_lvl[1]));

    std::cout << "testing_filter_entries"<< std::endl;
}
*/

inline void check_LU(Ilu_k_matlab_order& expected_precon, Iluk_moderncpp_standard_order & actual_precon) {
    for (int i = 0; i< expected_precon.L.size(); i++) {
        debug_assert(expected_precon.L[i], actual_precon.L[i].value,
            "expected " + std::to_string(expected_precon.L[i]) + ", got: " + std::to_string(actual_precon.L[i].value)
            + " at "+ std::to_string(actual_precon.L[i].row) + " "+ std::to_string(actual_precon.L[i].col));
    }

    for (int i = 0; i< expected_precon.U.size(); i++) {
        debug_assert(expected_precon.U[i], actual_precon.U[i].value,
            "expected " + std::to_string(expected_precon.U[i]) + ", got: " + std::to_string(actual_precon.U[i].value)
            + " at "+ std::to_string(actual_precon.L[i].row) + " "+ std::to_string(actual_precon.L[i].col));
    }

}

inline void correct_lvl0_decomposition_test() {

    auto [matrix, rhs, startvector]  = loadStandardTestingProblem();
    CSRMatrix A(matrix);
    std::vector<MatrixEntry> internal_LU = A.iluk(0);
    std::vector<MatrixEntry> test_A = matrix.getFullMatrix(false);

    auto expected_precon =  Ilu_k_matlab_order(matrix, 0);
    auto actual_precon = Iluk_moderncpp_standard_order(test_A,test_A.size(), 0);


    check_LU(expected_precon, actual_precon);
    std::cout << "passed correct_lvl0_decomposition_test"<< std::endl;
}

inline void correct_lvl1_decomposition_test() {

    auto [matrix, rhs, startvector]  = loadStandardTestingProblem();
    CSRMatrix A(matrix);
    std::vector<MatrixEntry> internal_LU = A.iluk(0);
    std::vector<MatrixEntry> test_A = matrix.getFullMatrix(false);

    auto expected_precon =  Ilu_k_matlab_order(matrix, 1);
    auto actual_precon = Iluk_moderncpp_standard_order(test_A,test_A.size(), 1);

    check_LU(expected_precon, actual_precon);
    std::cout << "passed correct_lvl1_decomposition_test"<< std::endl;
}

inline void correct_lvl2_decomposition_test() {

    auto [matrix, rhs, startvector]  = loadStandardTestingProblem();
    CSRMatrix A(matrix);
    std::vector<MatrixEntry> internal_LU = A.iluk(0);
    std::vector<MatrixEntry> test_A = matrix.getFullMatrix(false);

    auto expected_precon =  Ilu_k_matlab_order(matrix, 2);
    auto actual_precon = Iluk_moderncpp_standard_order(test_A,test_A.size(), 2);

    check_LU(expected_precon, actual_precon);
    std::cout << "passed correct_lvl2_decomposition_test"<< std::endl;
}

inline void correct_lvl3_decomposition_test() {

    auto [matrix, rhs, startvector]  = loadStandardTestingProblem();
    CSRMatrix A(matrix);
    std::vector<MatrixEntry> internal_LU = A.iluk(0);
    std::vector<MatrixEntry> test_A = matrix.getFullMatrix(false);

    auto expected_precon =  Ilu_k_matlab_order(matrix, 3);
    auto actual_precon = Iluk_moderncpp_standard_order(test_A,test_A.size(), 3);

    check_LU(expected_precon, actual_precon);
    std::cout << "passed correct_lvl3_decomposition_test"<< std::endl;
}




inline void iluk_moderncpp_standard_order_test() {
    //testing_grounds();
    //testing_sparse_vector_subtraction();//todo change to modified part
    //testing_filter_entries(); //todo change to modified part
    correct_lvl0_decomposition_test();
    correct_lvl1_decomposition_test();
    correct_lvl2_decomposition_test();
    correct_lvl3_decomposition_test();
}
#endif //ILUK_MODERNCPP_STANDARD_ORDER_TEST_H
