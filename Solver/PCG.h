//
// Created by jonas on 13.11.2024.
//

#ifndef PCG_H
#define PCG_H
#include <memory>

#include "..//Preconditioners/Preconditioner.h"
#include "../Matrix Management/DenseVector.h"
#include "../Matrix Management/COOMatrixSorted.h"

struct PcgSolution {
    int iterations;
    DenseVector solution_vector;
};

class PCG {

    public:
    PcgSolution calculate(const COOMatrixSorted& A,
                           DenseVector b,
                           const std::optional<DenseVector> &x_0 = std::nullopt,
                           const std::optional<std::shared_ptr<Preconditioner>> &precon = std::nullopt,
                           double tol = 1e-6,
                           int max_iter = 1000);

};



#endif //PCG_H
