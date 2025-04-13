//
// Created by jonas on 27.11.2024.
//

#ifndef COOMatrixSorted_H
#define COOMatrixSorted_H

#include <iostream>
#include <vector>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>

#include "CSCMatrix.h"
#include "CSRMatrix.h"
#include "DenseVector.h"
#include "SparseMatrix.h"

class CSRMatrix;
class CSCMatrix;

class COOMatrixSorted: public SparseMatrix {
private:
    std::vector<MatrixEntry> matrix;
    std::vector<std::pair<int, double*>> diags;

    void readFromMTX(const std::string& filepath);
    void get_diagonal_pointer_pairs();
    double binary_search(int target_row, int target_col);

public:
    int rows{}, cols{};

    // Constructor that reads a matrix from a .mtx file
    explicit COOMatrixSorted(const std::string& filepath);
    explicit COOMatrixSorted(const CSRMatrix& matrix);
    explicit COOMatrixSorted(const CSCMatrix& matrix);
    COOMatrixSorted(const std::vector<MatrixEntry>& matrix, int rows, int cols, bool column_sorted);

    DenseVector diagonal() override;
    double coeffRef(int row, int col) override;
    [[nodiscard]] std::vector<MatrixEntry> getStrictUpperTriangularPart(bool column_sorted) const override;
    [[nodiscard]] std::vector<MatrixEntry> getFullMatrix(bool column_sorted) const override;
    [[nodiscard]] std::vector<MatrixEntry> getStrictLowerTriangularPart(bool column_sorted) const override;

    // Matrix-vector multiplication: y = A * x
    DenseVector operator*(const DenseVector& x) const override ;
};

struct MatrixTriple {
    COOMatrixSorted matrix;
    DenseVector rhs;
    DenseVector startvector;

    explicit MatrixTriple(const std::string& path): matrix(path) {
    }
};
#endif //COOMatrixSorted_H
