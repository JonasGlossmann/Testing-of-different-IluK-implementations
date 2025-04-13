//
// Created by jonas on 03.04.2025.
//

#ifndef ILUK_BASICCPP_STANDARD_ORDER_H
#define ILUK_BASICCPP_STANDARD_ORDER_H
#include "Preconditioner.h"
#include "../Matrix Management/COOMatrixSorted.h"


class iluk_basiccpp_standard_order final :public Preconditioner{

    std::vector<double> L;
    std::vector<int> L_col;
    std::vector<int> L_row_ptr;
    std::vector<double> U;
    std::vector<int> U_col;
    std::vector<int> U_row_ptr;
    iluk_basiccpp_standard_order(COOMatrixSorted A_ext, int level);
    DenseVector apply(DenseVector residual) override;

};



#endif //ILUK_BASICCPP_STANDARD_ORDER_H
