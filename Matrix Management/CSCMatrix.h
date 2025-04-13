//
// Created by jonas on 08.12.2024.
//

#ifndef CSCMATRIX_H
#define CSCMATRIX_H
#include "COOMatrixSorted.h"
#include "CSRMatrix.h"
#include "SparseMatrix.h"

//todo testing
class COOMatrixSorted;
class CSRMatrix;

class CSCMatrix: public SparseMatrix {
private:

    void fromMatrixEntries(const std::vector<MatrixEntry>& matrix, int cols);

    //todo
public:
    std::vector<double> val;
    std::vector<int> rowIndex;
    std::vector<int> colPtr;

    int rows{}, cols{};
    explicit CSCMatrix(const std::string& filepath);
    explicit CSCMatrix(const COOMatrixSorted& matrix);
    explicit CSCMatrix(const CSRMatrix& matrix);
    CSCMatrix(std::vector<MatrixEntry> matrix, int rows, int cols);

    DenseVector diagonal() override;
    double coeffRef(int row, int col) override;

    [[nodiscard]] std::vector<MatrixEntry> getFullMatrix(bool column_sorted) const override;
    [[nodiscard]] std::vector<MatrixEntry> getStrictUpperTriangularPart(bool column_sorted) const override;
    [[nodiscard]] std::vector<MatrixEntry> getStrictLowerTriangularPart(bool column_sorted) const override;

    DenseVector operator*(const DenseVector &x) const override;
};



#endif //CSCMATRIX_H
