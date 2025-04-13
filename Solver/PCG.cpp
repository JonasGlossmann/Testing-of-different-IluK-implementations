//
// Created by jonas on 13.11.2024.
//

#include "PCG.h"

#include "../Matrix Management/DenseVector.h"
#include "../misc/TimeControl.h"
#include "../Preconditioners/IdentityPreconditioner.h"

PcgSolution PCG::calculate(const COOMatrixSorted& A,
                           DenseVector b,
                           const std::optional<DenseVector> &x_0,
                           const std::optional<std::shared_ptr<Preconditioner>> &precon,
                           double tol,
                           const int max_iter ) {
    std::shared_ptr<Preconditioner> preconditioner ;


    if (precon) {
        preconditioner = *precon;
    } else {
        preconditioner = std::make_shared<IdentityPreconditioner>();
    }

    DenseVector x;
    if(x_0) {
        x = x_0.value();
    } else {
        x = DenseVector(A.rows);
    }

    int n = 0;
    tol = tol* b.norm();
    DenseVector r = b - (A * x) ;
    DenseVector p = preconditioner->apply(r);
    double rsold = r.dot(p);


    while(n<max_iter) {
        if(r.norm() < tol) {
            break;
        }
        DenseVector Ap = A*p;
        double alpha = rsold/p.dot(Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        DenseVector h = preconditioner->apply(r);
        double rsnew = r.dot(h);
        p = h + (rsnew / rsold) * p;
        rsold = rsnew;
        n++;
    }

    PcgSolution solution;
    solution.iterations = n;
    solution.solution_vector = x;
    return solution;
}





