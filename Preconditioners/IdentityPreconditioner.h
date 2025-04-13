//
// Created by jonas on 13.11.2024.
//

#ifndef IDENTITYPRECONDITIONER_H
#define IDENTITYPRECONDITIONER_H
#include "Preconditioner.h"


class IdentityPreconditioner final :public Preconditioner {

public:
    DenseVector apply(DenseVector residual) override;
};



#endif //IDENTITYPRECONDITIONER_H
