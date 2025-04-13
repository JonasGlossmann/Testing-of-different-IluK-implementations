//
// Created by jonas on 28.11.2024.
//

#ifndef TESTMAIN_H
#define TESTMAIN_H

#include "MatrixTest/COOMatrixSorted_Test.h"
#include "MatrixTest/CSRMatrix.h"
#include "IluKTest/iluk_basiccpp_reverse_order_test.h"
#include "IluKTest/iluk_basiccpp_standard_order_test.h"
#include "IluKTest/iluk_moderncpp_reverse_order_test.h"
#include "IluKTest/iluk_moderncpp_standard_order_test.h"
#include "IluKTest/iluk_light_basiccpp_test.h"
#include "IluKTest/iluk_light_moderncpp_test.h"
#include "misc/Preconditioner_test.h"

inline void testMain() {
    //Preconditioner_test();

    //COOMatrixSorted_Test();
    //CSRMatrixTest();

    iluk_basiccpp_reverse_order_test();
    iluk_basiccpp_standard_order_test();
    iluk_moderncpp_reverse_order_test();
    iluk_moderncpp_standard_order_test();
    iluk_light_basiccpp_test();
    iluk_light_moderncpp_test();

}

#endif //TESTMAIN_H
