//
// Created by jonas on 13.11.2024.
//

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H
#include "../Matrix Management/DenseVector.h"
#include "../misc/data_structs.h"


class Preconditioner {

public:
//size of matrices and residual have to be matched already, otherwise this function may fail
static DenseVector ForwardBackwardSubstitutionModern(const DenseVector &residual,
                                                    const std::vector<MatrixEntry>& L,
                                                    const std::vector<MatrixEntry>& U) {
    DenseVector result = residual;

    //forward substitution
    for (auto&[row, col, val] : L) {
        result.data[row] -= val * result.data[col];
    }

    //backward substitution
    for (int idx = U.size() - 1; idx > -1; --idx) {
        if(auto &[row, col, val] = U[idx]; row != col) {
            result.data[row] -= val * result.data[col];
        } else {
            result.data[row] /= val;

        }
    }

    return result;
}

//size of matrices and residual have to be matched already, otherwise this function may fail
static DenseVector ForwardBackwardSubstitutionBasic(   const DenseVector& residual,
                                                        const std::vector<double> &L,
                                                        const std::vector<int> &L_row_ptr,
                                                        const std::vector<int> &L_col,
                                                        const std::vector<double> &U,
                                                        const std::vector<int> &U_row_ptr,
                                                        const std::vector<int> &U_col) {
    std::vector<double> result = residual.data;

    //foward substitution
    int num_rows = L_row_ptr.size() -1;
    for (int row = 0; row < num_rows; ++row) {
        int next_row =L_row_ptr[row+1];

        for(int idx = L_row_ptr[row]; idx< next_row; ++idx) {
            result[row] -= L[idx]*result[L_col[idx]];

        }
    }

    //backward substitution
    num_rows = U_row_ptr.size() -2;
    for(int row = num_rows; row >-1 ; --row) {
        int prev_row = U_row_ptr[row]-1;

        for(int idx = U_row_ptr[row+1]-1; idx > prev_row; --idx) {
            int col = U_col[idx];

            if(row != col) {
                result[row] -= U[idx] * result[col];
            } else {
                result[row] = result[row]/U[idx];
            }
        }
    }

    DenseVector result2(result); //can be removed, if pcg uses raw vectors instead of Class structure
    return result2;
}

//public:
    virtual ~Preconditioner() = default;

    virtual DenseVector apply(DenseVector residual) = 0;

};



#endif //PRECONDITIONER_H
