//
// Created by jonas on 06.12.2024.
//

#include "Ilu_k.h"
#include "../misc/TimeControl.h"


Ilu_k::Ilu_k(COOMatrixSorted A_ext, int level) {
    CSRMatrix A(A_ext);

    std::vector<MatrixEntry> LU =A.iluk(level) ;


    for (const auto& entry : LU) {
        if (entry.row > entry.col) {
            // Strictly lower triangular part goes to L
            this->L.emplace_back(entry.row, entry.col, entry.value);
        } else {
            // Upper triangular part goes to U
            this->U.emplace_back(entry.row, entry.col, entry.value);
        }
    }
    std::cout<<"Number of entries in preconditioning matrix:"<<LU.size() <<std::endl;

    std::ranges::sort(this->L, [](const MatrixEntry& a, const MatrixEntry& b) {
   return (a.row < b.row) || (a.row == b.row && a.col < b.col);
   });

}

DenseVector Ilu_k::apply(DenseVector residual) {
    return ForwardBackwardSubstitutionModern(residual, this->L, this->U);
}


