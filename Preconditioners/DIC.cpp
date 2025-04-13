//
// Created by jonas on 18.11.2024.
//

#include "DIC.h"
#include <iostream>

#include "../Matrix Management/MatrixManagement.h"

DIC::DIC(COOMatrixSorted A_ext): A(A_ext) {
    this->D = A_ext.diagonal();
    this->U_A = A_ext.getStrictUpperTriangularPart(true);

    for(auto& [row, col, val] : this->U_A) {
        this->D.data[col] -= val * val * this->D.data[row];
    }

    this->D = 1.0 / this->D;
}


DenseVector DIC::apply(DenseVector residual) {
    //Forward substitution
    DenseVector result = residual * this-> D;

    for(auto& [row, col, val] : this->U_A) {
        result.data[row] -= val * result.data[col] * this->D.data[row];

    }
    //backward substitution
    for (int idx = this->U_A.size() - 1; idx > -1; --idx) {
        auto &[row, col, val] = this->U_A[idx];

        result.data[col] -= this->D.data[col] * val * result.data[row];
    }

    return result;
}

