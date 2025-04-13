#include <iostream>

#include "Matrix Management/MatrixManagement.h"
#include "misc/TimeControl.h"
#include "Preconditioners/DIC.h"
#include "Preconditioners/Ilu_k.h"
#include "Solver/PCG.h"
#include "Tests/TestMain.h"

void runMatrix(MatrixTriple& triple, bool blank, bool dic, bool Ilu_0, bool Ilu_1, bool Ilu_2, bool Ilu_3) {
    PCG Solver;

    auto& [matrix, rhs, startvector] = triple;

    if (blank) {
        std::cout <<"No preconditioner:"  << std::endl;
        auto time_test = TimeControl();
        auto [iter, sol] = Solver.calculate(matrix, rhs, startvector);

        time_test.checkout();
        std::cout <<"Number of iterations: " <<iter << std::endl;
        std::cout << std::endl;
    }

    if (dic) {
        std::cout <<"DIC preconditioning:"  << std::endl;
        auto time_test = TimeControl();

        auto [iter, sol] = Solver.calculate(matrix, rhs, startvector,
                                              std::optional<std::shared_ptr<Preconditioner> >(new DIC(matrix)));
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
    }

    if (Ilu_0) {
        std::cout <<"ILU preconditioning with level of Fill = 0:"  << std::endl;
        auto time_test = TimeControl();
        PcgSolution result = Solver.calculate(matrix, rhs, startvector,
                                      std::optional<std::shared_ptr<Preconditioner> >(new Ilu_k(matrix,0)));
        auto& [iter, sol] = result;
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
        std::cout <<"ILU matlab preconditioning with level of Fill = 0:"  << std::endl;
        time_test = TimeControl();
        result = Solver.calculate(matrix, rhs, startvector,
                                      std::optional<std::shared_ptr<Preconditioner> >(new Ilu_k_matlab_order(matrix,0)));
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
    }

    if (Ilu_1) {
        std::cout <<"ILU preconditioning with level of Fill = 1:"  << std::endl;
        auto time_test = TimeControl();
        PcgSolution result = Solver.calculate(matrix, rhs, startvector,
                                      std::optional<std::shared_ptr<Preconditioner> >(new Ilu_k(matrix,1)));
        auto& [iter, sol] = result;
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
        std::cout <<"ILU matlab preconditioning with level of Fill = 1:"  << std::endl;
        time_test = TimeControl();
        result = Solver.calculate(matrix, rhs, startvector,
                                      std::optional<std::shared_ptr<Preconditioner> >(new Ilu_k_matlab_order(matrix,1)));
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
    }

    if (Ilu_2) {
        std::cout <<"ILU preconditioning with level of Fill = 2:"  << std::endl;
        auto time_test = TimeControl();
        PcgSolution result = Solver.calculate(matrix, rhs, startvector,
                                      std::optional<std::shared_ptr<Preconditioner> >(new Ilu_k(matrix,2)));
        auto& [iter, sol] = result;
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
        std::cout <<"ILU matlab preconditioning with level of Fill = 2:"  << std::endl;
        time_test = TimeControl();
        result = Solver.calculate(matrix, rhs, startvector,
                                      std::optional<std::shared_ptr<Preconditioner> >(new Ilu_k_matlab_order(matrix,2)));
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
    }

    if (Ilu_3) {
        std::cout <<"ILU preconditioning with level of Fill = 3:"  << std::endl;
        auto time_test = TimeControl();
        PcgSolution result = Solver.calculate(matrix, rhs, startvector,
                                      std::optional<std::shared_ptr<Preconditioner> >(new Ilu_k(matrix,3)));
        auto& [iter, sol] = result;
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
        std::cout <<"ILU matlab preconditioning with level of Fill = 3:"  << std::endl;
        time_test = TimeControl();
        result = Solver.calculate(matrix, rhs, startvector,
                                      std::optional<std::shared_ptr<Preconditioner> >(new Ilu_k_matlab_order(matrix,3)));
        time_test.checkout();
        std::cout <<"Number of Iterations: " <<iter << std::endl;
        std::cout << std::endl;
    }

}

int main() {
    const std::wstring root_directory = L"C:\\Users\\jonas\\Desktop\\Matrizen";
    const std::wstring mtx_filename = L"p_local_A.mtx";
    const std::wstring b_vec_filename = L"p_rhs.mtx";
    const std::wstring x_startvec_filename = L"p_initial_guess.mtx";

    //todo
    testMain();
    std::cout << "building matrix library" << std::endl;
    MatrixManagement Matrices(root_directory,
                             mtx_filename,
                             b_vec_filename,
                             x_startvec_filename);

    while (Matrices.validIndex()) {
        std::cout << "reading matrix" << std::endl;
        auto triple = Matrices.getMatrix();
        std::cout << "solving matrix \n" << std::endl;
        runMatrix(triple, false, false, true, true, true, true);

        system("pause");

        Matrices.next();
    }


    return 0;
}
