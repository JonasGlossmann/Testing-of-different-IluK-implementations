//
// Created by jonas on 08.12.2024.
//

#ifndef CSRMATRIX_H
#define CSRMATRIX_H
#include "SparseMatrix.h"
#include <algorithm>

#include "CSCMatrix.h"
#include "COOMatrixSorted.h"
#include "../misc/data_structs.h"
class CSCMatrix;
class COOMatrixSorted;

class CSRMatrix: public SparseMatrix {
    int rowBuilder;

    void kloop(int allowed_level, int i, std::vector<MatrixEntry> a_i, std::vector<MatrixLevel> a_i_lvl,
               std::vector<double> &new_diagonals, std::vector<MatrixEntry> &new_a,
               std::vector<int>& new_a_row_index, std::vector<MatrixLevel> &new_a_lvl);

    void jloop(
        int allowed_level, int kp1, std::vector<MatrixEntry> &a_i,
        std::vector<MatrixEntry> &a_k, std::vector<MatrixLevel> &a_i_lvl, std::vector<MatrixLevel> &a_k_lvl,
        std::vector<double> &new_diagonals, int a_ik_lvl, double a_ik);
    void getMatrixByMatrixEntries(std::vector<MatrixEntry>, int rows, int cols);
public:
    //todo move back
    std::vector<double> val;
    std::vector<int> rowPtr;
    std::vector<int> colIndex;
    int rows{}, cols{};

    explicit CSRMatrix(const std::string& filepath);
    explicit CSRMatrix(const COOMatrixSorted& matrix);
    explicit CSRMatrix(const CSCMatrix& matrix);
     CSRMatrix(std::vector<double> vals, std::vector<int> rowPtr, std::vector<int> colIndex);

    void addRow(std::vector<MatrixEntry> rowEntries);

    std::vector<MatrixEntry> getFullMatrix();

    CSRMatrix(const std::vector<MatrixEntry>& matrix, int rows, int cols);
    CSRMatrix(int rows, int cols);

    DenseVector diagonal() override;

    std::vector<double> standarddiagonal() const;

    double coeffRef(int row, int col) override;

    [[nodiscard]] std::vector<MatrixEntry> getStrictUpperTriangularPart(bool column_sorted) const override;
    [[nodiscard]] std::vector<MatrixEntry> getStrictLowerTriangularPart(bool column_sorted) const override;
    [[nodiscard]] std::vector<MatrixEntry> getFullMatrix(bool column_sorted) const override;

    DenseVector operator*(const DenseVector &x) const override;

    [[nodiscard]] std::vector<MatrixEntry> getRow(int row) const ;


    std::vector<MatrixEntry> iluk(int level);
};



#endif //CSRMATRIX_H
