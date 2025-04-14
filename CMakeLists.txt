cmake_minimum_required(VERSION 3.29)
project(Preconditioner_Testing)

set(CMAKE_CXX_STANDARD 26)

add_executable(Preconditioner_Testing main.cpp
        "Matrix Management/MatrixManagement.cpp"
        "Matrix Management/MatrixManagement.h"
        Solver/PCG.cpp
        Solver/PCG.h
        Preconditioners/Preconditioner.h
        Preconditioners/IdentityPreconditioner.cpp
        Preconditioners/IdentityPreconditioner.h
        Preconditioners/DIC.cpp
        Preconditioners/DIC.h
        misc/TimeControl.cpp
        misc/TimeControl.h
        "Matrix Management/COOMatrixSorted.cpp"
        "Matrix Management/COOMatrixSorted.h"
        "Matrix Management/DenseVector.cpp"
        "Matrix Management/DenseVector.h"
        Tests/MatrixTest/COOMatrixSorted_Test.h
        Tests/TestMain.h
        Preconditioners/Ilu_k.cpp
        Preconditioners/Ilu_k.h
        "Matrix Management/SparseMatrix.h"
        "Matrix Management/CSRMatrix.cpp"
        "Matrix Management/CSRMatrix.h"
        "Matrix Management/CSCMatrix.cpp"
        "Matrix Management/CSCMatrix.h"
        "Matrix Management/LduMatrix.cpp"
        "Matrix Management/LduMatrix.h"
        Tests/MatrixTest/CSRMatrix.h
        Preconditioners/Ilu_k_matlab_order.cpp
        Preconditioners/Ilu_k_matlab_order.h
        Tests/TestingTools/TestingTools.h
        Tests/IluKTest/iluk_moderncpp_standard_order_test.h
        Tests/IluKTest/iluk_moderncpp_reverse_order_test.h
        Tests/IluKTest/iluk_basiccpp_reverse_order_test.h
        Tests/IluKTest/iluk_basiccpp_standard_order_test.h
        Tests/IluKTest/iluk_light_basiccpp_test.h
        Tests/IluKTest/iluk_light_moderncpp_test.h
        Preconditioners/iluk_basiccpp_standard_order.cpp
        Preconditioners/iluk_basiccpp_standard_order.h
        Tests/misc/Preconditioner_test.h
        misc/data_structs.h
        Preconditioners/Iluk_moderncpp_standard_order.cpp
        Preconditioners/Iluk_moderncpp_standard_order.h
        Preconditioners/Ilu_k_moderncpp_alternative_standard_order.cpp
        Preconditioners/Ilu_k_moderncpp_alternative_standard_order.h
)
