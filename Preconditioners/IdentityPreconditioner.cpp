//
// Created by jonas on 13.11.2024.
//

#include "IdentityPreconditioner.h"

DenseVector IdentityPreconditioner::apply(DenseVector residual) {
    return residual;
}
