//
// Created by jonas on 08.12.2024.
//

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include "DenseVector.h"
#include "../misc/data_structs.h"

class SparseMatrix {
public:
    virtual ~SparseMatrix() = default;

private:
    virtual DenseVector diagonal() = 0;
    virtual double coeffRef(int row, int col) = 0;
    [[nodiscard]] virtual std::vector<MatrixEntry> getStrictUpperTriangularPart(bool column_sorted) const = 0;
    [[nodiscard]] virtual std::vector<MatrixEntry> getStrictLowerTriangularPart(bool column_sorted) const = 0;
    [[nodiscard]] virtual std::vector<MatrixEntry> getFullMatrix(bool column_sorted) const = 0;

    // Matrix-vector multiplication: y = A * x
    virtual DenseVector operator*(const DenseVector& x) const = 0;
};



#endif //SPARSEMATRIX_H
