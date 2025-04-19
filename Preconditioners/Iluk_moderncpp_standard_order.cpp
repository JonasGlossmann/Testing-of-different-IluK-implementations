//
// Created by jonas on 08.04.2025.
//

#include "Iluk_moderncpp_standard_order.h"

void copy_first_row(const std::vector<MatrixEntry>& A_old, std::vector<MatrixEntry>& A_new, std::vector<int>& A_lvl_new,
                    std::vector<size_t>& A_new_row_ptr, std::vector<MatrixEntry>& U, std::vector<double>& diags ) {
    for (size_t i = 0; i < A_old.size()&& A_old[i].row == 0; i++) {
        auto entry = A_old[i];
        A_new.push_back(entry);
        A_lvl_new.push_back(0);
        U.push_back(entry);
        if(entry.row == entry.col) {
            diags[0] = entry.value;
        }
    }
    A_new_row_ptr.push_back(0);
    A_new_row_ptr.push_back(A_new.size());
}

void get_row_by_index(int row_number, size_t& index, std::vector<MatrixEntry>&matrix, std::vector<int>& lvl_mtx,
                        std::vector<MatrixEntry>& row_return, std::vector<int>& row_lvl_return,
                        std::vector<double>& diags) {
    while (matrix[index].row==row_number) {
        auto entry = matrix[index];
        row_return.push_back(entry);
        row_lvl_return.push_back(lvl_mtx[index]);
        if(entry.row == entry.col) {
            diags[0] = entry.value;
        }
        ++index;
    }
}

void sparse_row_subtraction_with_level_and_starting_col (std::vector<MatrixEntry>& x,const std::vector<MatrixEntry>& y,
    std::vector<int>& x_lvl, const std::vector<int>& y_lvl, const double alpha,const int alpha_lvl,
     int& y_starting_index, int starting_col) {
    while(y_starting_index < y.size() &&y[y_starting_index].col <=starting_col) {
        ++y_starting_index;
    }
    int index = 0;
    for (int i = y_starting_index; i < y.size(); i++) {
        auto [y_row, y_col, y_value] = y[i];
        while (x[index].col < y_col) {
            ++index;
        }
        auto [x_row, x_col,x_value] = x[index];
        if(x_col == y_col) {
            x[index].value -= alpha*y_value;
            x_lvl[index] = std::min(x_lvl[index],alpha_lvl+y_lvl[i] + 1);
        }
        else {
            x.insert(x.begin()+index,MatrixEntry{x_row,y_col,-alpha*y_value});
            x_lvl.insert(x_lvl.begin()+index,alpha_lvl+y_lvl[i] + 1);
        }
    }
}

//@param A Matrix (need to be sorted by row)
//@param a_size number of cols/rows (matrix must be symmetrical)
//@return
Iluk_moderncpp_standard_order::Iluk_moderncpp_standard_order(std::vector<MatrixEntry>& A, size_t a_size, const int level) {

    std::vector<int> A_lvl_old(A.size(),0);
    std::vector<MatrixEntry> A_new;
    std::vector<int> A_lvl_new;
    std::vector<size_t> A_new_row_ptr;
    std::vector<double> diags(a_size);

    copy_first_row(A,A_new,A_lvl_new,A_new_row_ptr, this->U, diags);

    size_t index = 0;
    //iloop
    for (int i = 1; i < a_size; ++i) {
        std::vector<MatrixEntry> a_i;
        std::vector<int> a_i_lvl;

        get_row_by_index(i,index,A,A_lvl_old,a_i,a_i_lvl,diags);

        size_t vec_index = 0;
        int last_col = 0;
        int y_index = 0;

/*
        //kloop
        while( a_i[vec_index].col<i && vec_index<a_i.size()) {
            auto& [row,col,value] = a_i[vec_index];

            if(col>last_col &&a_i_lvl[vec_index]<=level) {
                last_col = col;
                double alpha = value/diags[last_col];
                this->L.push_back({row,col,alpha});

                std::vector<MatrixEntry> a_k(A_new.begin()+A_new_row_ptr[col],
                    A_new.begin()+A_new_row_ptr[col]);
                std::vector<int> a_k_lvl(A_lvl_new.begin() + A_new_row_ptr[col],
                    A_lvl_new.begin()+A_new_row_ptr[col]);

                sparse_row_subtraction_with_level_and_starting_col(a_i,a_k ,a_i_lvl,a_k_lvl,
                    alpha,a_i_lvl[vec_index],y_index,col);
            }
            vec_index++;
        }

*/
    }


}