//
// Created by jonas on 05.04.2025.
//

#ifndef PRECONDITIONERTESTING_H
#define PRECONDITIONERTESTING_H

//Testing Forward and backward Substitution against examplatory implementation
inline void testing_forward_backward_substitution_basic() {

    auto [matrix, rhs, startvector] = loadStandardTestingProblem();

    Ilu_k_matlab_order control_preconditioner(matrix, 0);
    DenseVector control_result = control_preconditioner.apply(rhs);

    std::vector<double> L = control_preconditioner.L;
    std::vector<double> U = control_preconditioner.U;
    std::vector<int> L_row_ptr = control_preconditioner.L_row_ptr;
    std::vector<int> U_row_ptr = control_preconditioner.U_row_ptr;
    std::vector<int> L_col = control_preconditioner.L_col;
    std::vector<int> U_col = control_preconditioner.U_col;

    DenseVector result = Preconditioner::ForwardBackwardSubstitutionBasic(rhs,L,L_row_ptr,L_col,U,U_row_ptr,U_col);

    for (int i=0;i<control_result.size;i++) {
        assert(control_result.data[i]==result.data[i]);
    }
    std::cout<<"testing_forward_backward_substitution_basic passed"<<std::endl;
};

inline void testing_forward_backward_substitution_modern() {
    auto [matrix, rhs, startvector] = loadStandardTestingProblem();
    Ilu_k_matlab_order control_preconditioner(matrix, 0);
    DenseVector control_result = control_preconditioner.apply(rhs);

    Ilu_k preconditioner(matrix, 0);

    DenseVector result = Preconditioner::ForwardBackwardSubstitutionModern(rhs,preconditioner.L,preconditioner.U);

    for (int i=0;i<control_result.size;i++) {
        assert(control_result.data[i]==result.data[i]);
    }
    std::cout << "testing_forward_backward_substitution_modern passed" << std::endl;
}

inline void Preconditioner_test() {
    testing_forward_backward_substitution_basic();
    testing_forward_backward_substitution_modern();

};
#endif //PRECONDITIONERTESTING_H
