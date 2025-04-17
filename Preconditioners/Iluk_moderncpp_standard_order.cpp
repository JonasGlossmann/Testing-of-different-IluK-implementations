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


    std::vector<int> lvl_old(A.size(),0);
    std::vector<int> lvl_new;
    std::vector<MatrixEntry> new_entries;
    std::vector<double> diag(a_size);
    std::vector<int> new_entries_row_indices(a_size+1,0);

    size_t index = 0;
    while (A[index].row == 0) {
        lvl_new.push_back(lvl_old[index]);
        new_entries.push_back(A[index]);
        if(A[index].col == 0) {
            diag[0] = A[index].value;
        }
        ++index;
    }
    new_entries_row_indices[1] = new_entries.size();

    for (size_t i = 1; i < a_size; i++) {
        std::vector<MatrixEntry> a_i ;
        std::vector<int> a_i_lvl;

        while (A[index].row == i) {
            a_i.push_back(A[index]);
            a_i_lvl.push_back(lvl_old[index]);
            if(A[index].col == i) {
                diag[i] = A[index].value;
            }
            ++index;
        }

        size_t vector_index = 0;
        size_t vector_col = 0;

        while (a_i[vector_index].row == i) {
            auto& [a_i_row, a_i_col, a_i_value] = a_i[vector_index];
            if(a_i_col >=vector_col && a_i_col < i && a_i_lvl[vector_index] <= level ) {

                a_i_value = a_i_value/diag[a_i_col];
                double alpha = a_i_value;
                int alpha_lvl = a_i_lvl[vector_index];

                std::vector<MatrixEntry> a_k(new_entries.begin()+new_entries_row_indices[a_i_col],
                    new_entries.begin()+new_entries_row_indices[a_i_col]);
                std::vector<int> a_k_lvl(lvl_new.begin() + new_entries_row_indices[a_i_col],
                    lvl_new.begin()+new_entries_row_indices[a_i_col]);

                sparse_row_subtraction_with_level(a_i,a_k ,a_i_lvl,a_k_lvl ,alpha,alpha_lvl);

            }
            ++vector_index;
        }


        filter_entries(a_i,a_i_lvl,level);
        new_entries.insert(new_entries.end(),a_i.begin(),a_i.end());
        lvl_new.insert(lvl_new.end(),a_i_lvl.begin(),a_i_lvl.end());
        new_entries_row_indices[i+1] = new_entries.size();
    }

    //separate L/U matrix into L and U
    for (const auto& entry : new_entries) {
        if (entry.row > entry.col) {
            // Strictly lower triangular part goes to L
            this->L.emplace_back(entry.row, entry.col, entry.value);
        } else {
            // Upper triangular part goes to U
            this->U.emplace_back(entry.row, entry.col, entry.value);
        }
    }


}