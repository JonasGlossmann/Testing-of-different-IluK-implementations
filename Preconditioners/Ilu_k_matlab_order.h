//
// Created by jonas on 11.03.2025.
//

#ifndef ILU_K_MATLAB_ORDER_H
#define ILU_K_MATLAB_ORDER_H
#include "Preconditioner.h"
#include "../Matrix Management/COOMatrixSorted.h"

class Ilu_k_matlab_order final :public Preconditioner{
private:


public:

    //todo move back
    std::vector<double> L;
    std::vector<int> L_col;
    std::vector<int> L_row_ptr;
    std::vector<double> U;
    std::vector<int> U_col;
    std::vector<int> U_row_ptr;
    Ilu_k_matlab_order(COOMatrixSorted A_ext, int level);
    DenseVector apply(DenseVector residual) override;
};




#endif //ILU_K_MATLAB_ORDER_H
