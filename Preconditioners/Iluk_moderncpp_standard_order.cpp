//
// Created by jonas on 08.04.2025.
//

#include "Iluk_moderncpp_standard_order.h"

void copy_first_row(const std::vector<MatrixEntry>& A_old, std::vector<MatrixEntry>& A_new, std::vector<int>& A_lvl_new,
                    std::vector<size_t>& A_new_row_ptr, std::vector<MatrixEntry>& U ) {
    for (size_t i = 0; i < A_old.size()&& A_old[i].row == 0; i++) {
        A_new.push_back(A_old[i]);
        A_lvl_new.push_back(0);
        U.push_back(A_old[i]);
    }
    A_new_row_ptr.push_back(0);
    A_new_row_ptr.push_back(A_new.size());
}

void get_row_by_index(int row_number, size_t& index, std::vector<MatrixEntry>&matrix, std::vector<int>& lvl_mtx,
                        std::vector<MatrixEntry>& row_return, std::vector<int>& row_lvl_return) {
    std::vector<MatrixEntry> row;
    while (matrix[index].row==row_number) {
        row_return.push_back(matrix[index]);
        row_lvl_return.push_back(lvl_mtx[index]);
        ++index;
    }
}

//@param A Matrix (need to be sorted by row)
//@param a_size number of cols/rows (matrix must be symmetrical)
//@return
Iluk_moderncpp_standard_order::Iluk_moderncpp_standard_order(std::vector<MatrixEntry>& A, size_t a_size, int level) {

    std::vector<int> A_lvl_old(A.size(),0);
    std::vector<MatrixEntry> A_new;
    std::vector<int> A_lvl_new;
    std::vector<size_t> A_new_row_ptr;

    copy_first_row(A,A_new,A_lvl_new,A_new_row_ptr, this->U);

    size_t index = 0;
    for (int i = 1; i<a_size; ++i) {
        std::vector<MatrixEntry> a_i;
        std::vector<int> a_i_lvl;

        get_row_by_index(i,index,A,A_lvl_old,a_i,a_i_lvl);

        size_t vec_index = 0;
        int last_col = 0;
        while(a_i[vec_index].col<i &&a_i_lvl[vec_index]<=level) {

            if(a_i[vec_index].col>last_col) {
                last_col = a_i[vec_index].col;
                double alpha =
            }
        }


    }


}