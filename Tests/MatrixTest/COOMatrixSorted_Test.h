//
// Created by jonas on 28.11.2024.
//

#ifndef COOMatrixSorted_TEST_H
#define COOMatrixSorted_TEST_H

#include "../../Matrix Management/COOMatrixSorted.h"
#include "../../Matrix Management/DenseVector.h"
#include <cassert>


inline COOMatrixSorted getMatrix() {
    return COOMatrixSorted(R"(..\Tests\TestFiles\testfile.mtx)");
}

inline DenseVector getVector() {
    DenseVector testvector(3);
    testvector.data[0] = 1;
    testvector.data[1] = 2;
    testvector.data[2] = 3;
    return testvector;
}

inline void constructorTest() {
    COOMatrixSorted matrix = getMatrix();

    assert(matrix.rows == 3);
    assert(matrix.cols == 3);

    assert(matrix.coeffRef(0, 0) == 1);
    assert(matrix.coeffRef(0, 1) == 2);
    assert(matrix.coeffRef(0, 2) == 3);
    assert(matrix.coeffRef(1, 0) == 4);
    assert(matrix.coeffRef(1, 1) == 5);
    assert(matrix.coeffRef(1, 2) == 6);
    assert(matrix.coeffRef(2, 0) == 0);
    assert(matrix.coeffRef(2, 1) == 7);
    assert(matrix.coeffRef(2, 2) == 8);
    std::cout << "COOMatrixSorted_Test constructor test passed" << std::endl;
}

inline void MatrixVectorMultiplicationTest_1() {
    COOMatrixSorted matrix = getMatrix();
    DenseVector testvector(3);
    testvector.data[0] = 1;
    auto result = matrix * testvector;
    assert(result.data[0] == 1);
    assert(result.data[1] == 4);
    assert(result.data[2] == 0);

    std::cout << "COOMatrixSorted_Test MatrixVectorMultiplication_1 test passed" << std::endl;
}

inline void MatrixVectorMultiplicationTest_2() {
    COOMatrixSorted matrix = getMatrix();
    DenseVector testvector(3);
    testvector.data[1] = 1;
    auto result = matrix * testvector;
    assert(result.data[0] == 2);
    assert(result.data[1] == 5);
    assert(result.data[2] == 7);

    std::cout << "COOMatrixSorted_Test MatrixVectorMultiplication_2 test passed" << std::endl;
}

inline void MatrixVectorMultiplicationTest_3() {
    COOMatrixSorted matrix = getMatrix();
    DenseVector testvector(3);
    testvector.data[2] = 1;
    auto result = matrix * testvector;
    assert(result.data[0] == 3);
    assert(result.data[1] == 6);
    assert(result.data[2] == 8);

    std::cout << "COOMatrixSorted_Test MatrixVectorMultiplication_3 test passed" << std::endl;
}

inline void getFullMatrixTest_1() {
    COOMatrixSorted matrix = getMatrix();
    auto result = matrix.getFullMatrix(true);

    auto entry = result[0];
    assert(entry.row == 0);
    assert(entry.col == 0);
    assert(entry.value == 1.0);
    entry = result[1];
    assert(entry.row == 1);
    assert(entry.col == 0);
    assert(entry.value == 4.0);
}

inline void getFullMatrixTest_2() {
    COOMatrixSorted matrix = getMatrix();
    auto result = matrix.getFullMatrix(false);

    auto entry = result[0];
    assert(entry.row == 0);
    assert(entry.col == 0);
    assert(entry.value == 1.0);
    entry = result[1];
    assert(entry.row == 0);
    assert(entry.col == 1);
    assert(entry.value == 2.0);
}

inline void COOMatrixSorted_Test() {
    constructorTest();
    MatrixVectorMultiplicationTest_1();
    MatrixVectorMultiplicationTest_2();
    MatrixVectorMultiplicationTest_3();
    getFullMatrixTest_1();
    getFullMatrixTest_2();
}

#endif //COOMatrixSorted_TEST_H
