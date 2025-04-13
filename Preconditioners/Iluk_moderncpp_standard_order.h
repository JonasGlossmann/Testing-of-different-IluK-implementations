//
// Created by jonas on 08.04.2025.
//

#ifndef ILUK_MODERNCPP_STANDARD_ORDER_H
#define ILUK_MODERNCPP_STANDARD_ORDER_H
#include "Preconditioner.h"
#include <iostream>
#include "../Matrix Management/COOMatrixSorted.h"


void sparse_row_subtraction_with_level (std::vector<MatrixEntry>& x,const std::vector<MatrixEntry>& y,
     std::vector<int>& x_lvl, const std::vector<int>& y_lvl, double alpha, int alpha_lvl);

void filter_entries(std::vector<MatrixEntry>& vector,std::vector<int>& lvl, int max_lvl);


class Iluk_moderncpp_standard_order final: public Preconditioner {
public:
    std::vector<MatrixEntry> L;
    std::vector<MatrixEntry> U;

    Iluk_moderncpp_standard_order(COOMatrixSorted A_ext, int level);
    DenseVector apply(DenseVector residual) override;
};



#endif //ILUK_MODERNCPP_STANDARD_ORDER_H
