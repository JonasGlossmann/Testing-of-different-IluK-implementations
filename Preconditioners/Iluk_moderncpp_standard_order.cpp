//
// Created by jonas on 08.04.2025.
//

#include "Iluk_moderncpp_standard_order.h"

void copy_first_row(const std::vector<MatrixEntry>& A_old, std::vector<MatrixEntry>& A_new, std::vector<int>& A_lvl_new,
                    std::vector<size_t>& A_new_row_ptr, std::vector<MatrixEntry>& U, std::vector<double>& diags,
                    size_t& index) {
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
    index = A_new_row_ptr[1];
}

void get_row_by_index(int row_number, size_t& index, std::vector<MatrixEntry>&matrix, std::vector<int>& lvl_mtx,
                        std::vector<MatrixEntry>& row_return, std::vector<int>& row_lvl_return) {
    while (matrix[index].row==row_number) {
        auto entry = matrix[index];
        row_return.push_back(entry);
        row_lvl_return.push_back(lvl_mtx[index]);
        ++index;
    }
}


void sparse_row_subtraction_with_level_and_starting_col (std::vector<MatrixEntry>& a_i,const std::vector<MatrixEntry>& a_k,
    std::vector<int>& a_i_lvl, const std::vector<int>& a_k_lvl, const double alpha,const int alpha_lvl,
     int& a_k_starting_index, int k, std::vector<double>& diags) {

    size_t y_size = a_k.size();

    while(a_k_starting_index < y_size &&a_k[a_k_starting_index].col <=k) {
        ++a_k_starting_index;
    }

    int index = 0;

    for (int i = a_k_starting_index; i < y_size; i++) {
        auto [y_row, y_col, y_value] = a_k[i];
        while (a_i[index].col < y_col) {
            ++index;
        }
        auto [x_row, x_col,x_value] = a_i[index];

        if(x_col == y_col) {
            a_i[index].value -= alpha*y_value;
            a_i_lvl[index] = std::min(a_i_lvl[index],alpha_lvl+a_k_lvl[i] + 1);
            if(x_row == x_col) {
                diags[x_row] =  a_i[index].value;
            }
            //todo debug
            if(x_row == 102 && x_col == 101) {
                std::cout << x_row << ' ' << x_col << ' ' << x_value << " " << alpha << " " << y_value << std::endl;
            }
            //todo
        }
        else {
            const auto tmp = -alpha*y_value;
            a_i.insert(a_i.begin()+index,MatrixEntry{x_row,y_col,tmp});
            a_i_lvl.insert(a_i_lvl.begin()+index,alpha_lvl+a_k_lvl[i] + 1);
            if(x_row == x_col) {
                diags[x_row] =  tmp;
            }
            //todo debug
            if(x_row == 102 && x_col == 101) {
                std::cout << x_row << ' ' << x_col << ' ' << x_value << " " << alpha << " " << y_value << std::endl;
            }
            //todo
        }
    }
}

//@param A Matrix (need to be sorted by row)
//@param a_size number of cols/rows (matrix must be symmetrical)
Iluk_moderncpp_standard_order::Iluk_moderncpp_standard_order(std::vector<MatrixEntry>& A,const size_t a_size, const int level) {

    std::vector<int> A_lvl_old(A.size(),0);
    std::vector<MatrixEntry> A_new;
    std::vector<int> A_lvl_new;
    std::vector<size_t> A_new_row_ptr;
    std::vector<double> diags(a_size);

    size_t obtain_a_i_index = 0; //index for obtaining a_i without the need to iterate through A all the time
    copy_first_row(A,A_new,A_lvl_new,A_new_row_ptr, this->U, diags, obtain_a_i_index);


    //iloop
    for (int i = 1; i < a_size; ++i) {
        std::vector<MatrixEntry> a_i; //row i of matrix A (old)
        std::vector<int> a_i_lvl; //lvl of row i
        std::vector<MatrixEntry> filtered_a_i; //row i of matrix A (old)
        std::vector<int> filtered_a_i_lvl; //lvl of row i

        //initializes a_i and a_i_lvl with row i of A (old)
        get_row_by_index(i,obtain_a_i_index,A,A_lvl_old,a_i,a_i_lvl);


        size_t a_i_index = 0; // current position in vector a_i
        int j_index = 0; //starting index for j-loop

        //kloop
        while( a_i[a_i_index].col<i && a_i_index<a_i.size()) {
            auto& [row,col,value] = a_i[a_i_index];

            if(a_i_lvl[a_i_index]<=level) {

                double alpha = value/diags[col];
                value = alpha;
                auto entry = MatrixEntry{row,col,value};
                this->L.push_back(entry);
                filtered_a_i.push_back(entry);
                filtered_a_i_lvl.push_back(a_i_lvl[a_i_index]);

                std::vector<MatrixEntry> a_k(A_new.begin()+A_new_row_ptr[col],
                    A_new.begin()+A_new_row_ptr[col+1]);
                std::vector<int> a_k_lvl(A_lvl_new.begin() + A_new_row_ptr[col],
                    A_lvl_new.begin()+A_new_row_ptr[col+1]);

                sparse_row_subtraction_with_level_and_starting_col(a_i,a_k ,a_i_lvl,a_k_lvl,
                    alpha,a_i_lvl[a_i_index],j_index,col, diags);
            }
            a_i_index++;
        }

        //filtering remaining entries of a_i and write them to U and to the new Matrix A_new (needed for remaining
        //iterations)
        while(a_i_index<a_i.size()) {
            if(a_i_lvl[a_i_index]<=level) {
                MatrixEntry entry = a_i[a_i_index];
                U.push_back(entry);
                filtered_a_i.push_back(entry);
                filtered_a_i_lvl.push_back(a_i_lvl[a_i_index]);
            }
            ++a_i_index;
        }

        A_new.insert(A_new.end(),filtered_a_i.begin(),filtered_a_i.end());
        A_lvl_new.insert(A_lvl_new.end(),filtered_a_i_lvl.begin(),filtered_a_i_lvl.end());
        A_new_row_ptr.push_back(A_new.size());
    }
}

DenseVector Iluk_moderncpp_standard_order::apply(DenseVector residual) {
    return Iluk_moderncpp_standard_order::ForwardBackwardSubstitutionModern(residual, this->L, this->U);
}
