//
// Created by jonas on 18.11.2024.
//

#ifndef DIC_H
#define DIC_H
#include "Preconditioner.h"
#include "../Matrix Management/COOMatrixSorted.h"

class DIC final:public Preconditioner {
    DenseVector D;
    COOMatrixSorted A;
    std::vector<MatrixEntry> U_A; //initially sorted ascending by columns

public:
    DIC(COOMatrixSorted A_ext);
    DenseVector apply(DenseVector residual) override;
};



#endif //DIC_H
