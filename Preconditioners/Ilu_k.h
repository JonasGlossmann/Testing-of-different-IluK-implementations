//
// Created by jonas on 06.12.2024.
//

#ifndef ILU_K_H
#define ILU_K_H
#include "Preconditioner.h"
#include "../Matrix Management/COOMatrixSorted.h"

class Ilu_k final :public Preconditioner {


public:
    std::vector<MatrixEntry> L;
    std::vector<MatrixEntry> U;
    Ilu_k(COOMatrixSorted A_ext, int level);
    DenseVector apply(DenseVector residual) override;
};



#endif //ILU_K_H
