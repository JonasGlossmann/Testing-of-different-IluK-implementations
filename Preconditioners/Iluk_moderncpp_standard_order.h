//
// Created by jonas on 08.04.2025.
//

#ifndef ILUK_MODERNCPP_STANDARD_ORDER_H
#define ILUK_MODERNCPP_STANDARD_ORDER_H
#include "Preconditioner.h"
#include "../misc/data_structs.h"
#include "../Matrix Management/COOMatrixSorted.h"



class Iluk_moderncpp_standard_order final: public Preconditioner {
public:
    std::vector<MatrixEntry> L;
    std::vector<MatrixEntry> U;


    Iluk_moderncpp_standard_order(std::vector<MatrixEntry>& A_ext, size_t a_cols, int level);
    DenseVector apply(DenseVector residual) override;
};



#endif //ILUK_MODERNCPP_STANDARD_ORDER_H
