//
// Created by jonas on 14.04.2025.
//

#include "Ilu_k_moderncpp_alternative_standard_order.h"

//calculates x_i = x_i - alpha * y_i
void sparse_row_subtraction_with_level (std::vector<MatrixEntry>& x,const std::vector<MatrixEntry>& y,
    std::vector<int>& x_lvl, const std::vector<int>& y_lvl, const double alpha,const int alpha_lvl) {
    unsigned int index = 0;
    for (int i = 0; i < y.size(); i++) {
        auto [y_row, y_col, y_value] = y[i];
        while (x[index].col < y_col) {
            ++index;
        }
        if(x[index].col == y_col) {
            x[index].value -= alpha*y_value;
            x_lvl[index] = std::min(x_lvl[index],alpha_lvl+y_lvl[i] + 1);
        }
        else {
            MatrixEntry new_entry = {x[index].row,y_col,-alpha*y_value};
            x.insert(x.begin()+index,new_entry);
            x_lvl.insert(x_lvl.begin()+index,alpha_lvl+y_lvl[i] + 1);
        }
    }
}

//takes a vector and the corresponding levels of each entry and deletes every vector and level entry
//with a level value over max_lvl
void filter_entries(std::vector<MatrixEntry>& vector, std::vector<int>& lvl,const int max_lvl) {
    std::vector<MatrixEntry> new_vector;
    std::vector<int> new_lvl;

    for (int i = 0; i < vector.size(); i++) {
        if(lvl[i]<=max_lvl) {
            new_vector.push_back(vector[i]);
            new_lvl.push_back(lvl[i]);
        }
    }
    vector = new_vector;
    lvl = new_lvl;
}


Ilu_k_moderncpp_alternative_standard_order::Ilu_k_moderncpp_alternative_standard_order(COOMatrixSorted A_ext, int level) {
    std::vector<MatrixEntry> old_entries = A_ext.getFullMatrix(false);
    std::vector<int> lvl_old(old_entries.size(),0);
    std::vector<int> lvl_new;
    std::vector<MatrixEntry> new_entries;
    std::vector<double> diag(A_ext.rows);
    std::vector<int> new_entries_row_indices(A_ext.rows+1,0);

    size_t index = 0;
    while (old_entries[index].row == 0) {
        lvl_new.push_back(lvl_old[index]);
        new_entries.push_back(old_entries[index]);
        if(old_entries[index].col == 0) {
            diag[0] = old_entries[index].value;
        }
        ++index;
    }
    new_entries_row_indices[1] = new_entries.size();

    for (size_t i = 1; i < A_ext.rows; i++) {
        std::vector<MatrixEntry> a_i ;
        std::vector<int> a_i_lvl;

        while (old_entries[index].row == i) {
            a_i.push_back(old_entries[index]);
            a_i_lvl.push_back(lvl_old[index]);
            if(old_entries[index].col == i) {
                diag[i] = old_entries[index].value;
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

DenseVector Ilu_k_moderncpp_alternative_standard_order::apply(DenseVector residual) {
    return ForwardBackwardSubstitutionModern(residual, this->L, this->U);
}