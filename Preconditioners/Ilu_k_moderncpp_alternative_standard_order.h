//
// Created by jonas on 14.04.2025.
//

#ifndef ILU_K_MODERNCPP_ALTERNATIVE_STANDARD_ORDER_H
#define ILU_K_MODERNCPP_ALTERNATIVE_STANDARD_ORDER_H

#include "Preconditioner.h"
#include "../misc/data_structs.h"
#include "../Matrix Management/COOMatrixSorted.h"

void sparse_row_subtraction_with_level (std::vector<MatrixEntry>& x,const std::vector<MatrixEntry>& y,
     std::vector<int>& x_lvl, const std::vector<int>& y_lvl, double alpha, int alpha_lvl);

void filter_entries(std::vector<MatrixEntry>& vector,std::vector<int>& lvl, int max_lvl);


class Ilu_k_moderncpp_alternative_standard_order final: public Preconditioner {
public:
    std::vector<MatrixEntry> L;
    std::vector<MatrixEntry> U;

    Ilu_k_moderncpp_alternative_standard_order(COOMatrixSorted A_ext, int level);
    DenseVector apply(DenseVector residual) override;
};


#endif //ILU_K_MODERNCPP_ALTERNATIVE_STANDARD_ORDER_H
