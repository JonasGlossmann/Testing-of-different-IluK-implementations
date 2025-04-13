//
// Created by jonas on 08.12.2024.
//

#include "CSRMatrix.h"

#include <ranges>
#include <utility>

#include "../misc/TimeControl.h"

void CSRMatrix::getMatrixByMatrixEntries(std::vector<MatrixEntry> matrix, const int rows, const int cols) {
    this->rows = rows;
    this->cols = cols;

    std::ranges::sort(matrix, [](const MatrixEntry& a, const MatrixEntry& b) {
    return (a.row < b.row) || (a.row == b.row && a.col < b.col);
    });

    this->rowPtr.resize(rows + 1, 0);
    std::vector<double> values;

    int current_row = -1;
    for (const auto& entry : matrix) {
        if (entry.row != current_row) {
            for (int i = current_row + 1; i <= entry.row; ++i) {
                this->rowPtr[i] = this->colIndex.size();
            }
            current_row = entry.row;
        }
        this->colIndex.push_back(entry.col);
        this->val.push_back(entry.value);
    }

    // Finalize row_ptr for remaining rows
    for (int i = current_row + 1; i <= rows; ++i) {
        this->rowPtr[i] = this->colIndex.size();
    }
}

CSRMatrix::CSRMatrix(const std::vector<MatrixEntry>& matrix, const int rows, const int cols) {
   getMatrixByMatrixEntries(matrix, rows, cols);
}

CSRMatrix::CSRMatrix(const int rows, const int cols) {
    this->rows = rows;
    this->cols = cols;
    this->rowBuilder = 0;
    this->rowPtr.resize(rows + 1, 0);
}

CSRMatrix::CSRMatrix(std::vector<double> vals, std::vector<int> rowPtr, std::vector<int> colIndex) {
    this->rowPtr = std::move(rowPtr);
    this->val = std::move(vals);
    this->colIndex = std::move(colIndex);
    this->rows = this->rowPtr.size()-1;
    this->cols = this->rowPtr.size()-1;
    this->rowBuilder = 0;
}

CSRMatrix::CSRMatrix(const COOMatrixSorted& matrix) {
    std::vector<MatrixEntry> matrixEntries = matrix.getFullMatrix(false);
    getMatrixByMatrixEntries(matrixEntries, matrix.rows, matrix.cols);
}

CSRMatrix::CSRMatrix(const std::string &filepath) {
    COOMatrixSorted matrix(filepath);
    std::vector<MatrixEntry> matrixEntries = matrix.getFullMatrix(false);
    getMatrixByMatrixEntries(matrixEntries, matrix.rows, matrix.cols);
}

CSRMatrix::CSRMatrix(const CSCMatrix &matrix) {
    COOMatrixSorted mat(matrix);
    std::vector<MatrixEntry> matrixEntries = matrix.getFullMatrix(false);
    getMatrixByMatrixEntries(matrixEntries, matrix.rows, matrix.cols);
}


void CSRMatrix::addRow(std::vector<MatrixEntry> rowEntries) {
    if(rowEntries[0].row != this->rowBuilder) {
        throw std::invalid_argument("Matrix row does not match necessary row");
    }

}

DenseVector CSRMatrix::diagonal() {
    if (rows != cols) {
        throw std::invalid_argument("Matrix must be square to extract diagonal.");
    }

    DenseVector diagonal(rows); // Initialize diagonal with zeros

    for (int row = 0; row < rows; ++row) {
        for (int idx = this->rowPtr[row]; idx < this->rowPtr[row + 1]; ++idx) {
            if (this->colIndex[idx] == row) { // Check if this is a diagonal element
                diagonal.data[row] = this->val[idx];
                break; // No need to continue searching this row
            }
        }
    }

    return diagonal;
}

std::vector<double> CSRMatrix::standarddiagonal() const {
    std::vector<double> diagonal(rows);

    for (int row = 0; row < rows; ++row) {
        for (int idx = this->rowPtr[row]; idx < this->rowPtr[row + 1]; ++idx) {
            if (this->colIndex[idx] == row) { // Check if this is a diagonal element
                diagonal[row] = this->val[idx];
                break; // No need to continue searching this row
            }
        }
    }
    return diagonal;
}

double CSRMatrix::coeffRef(const int row, const int col) {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("Row or column index out of bounds.");
    }


    // Search for column `j` in the current row
    for (int idx = this->rowPtr[row]; idx < this->rowPtr[row + 1]; ++idx) {
        if (this->colIndex[idx] == col) {
            return this->val[idx]; // Found the element
        }
    }

    return 0.0; // Element not found
}

std::vector<MatrixEntry> CSRMatrix::getFullMatrix(const bool column_sorted) const {
    std::vector<MatrixEntry> entries;

    // Extract all non-zero entries
    for (int row = 0; row < rows; ++row) {
        for (int idx = this->rowPtr[row]; idx < this->rowPtr[row + 1]; ++idx) {
            entries.emplace_back(row, this->colIndex[idx], this->val[idx]);
        }
    }

    // Sort by column if requested
    if (column_sorted) {
        std::ranges::sort(entries, [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.col < b.col) || (a.col == b.col && a.row < b.row);
        });
    }

    return entries;
}

std::vector<MatrixEntry> CSRMatrix::getStrictLowerTriangularPart(const bool column_sorted) const {
    std::vector<MatrixEntry> entries;

    for (int row = 0; row < rows; ++row) {
        for (int idx = this->rowPtr[row]; idx < this->rowPtr[row + 1]; ++idx) {
            if (int col = this->colIndex[idx]; col < row) { // Strict lower triangular condition
                entries.emplace_back(row, col, this->val[idx]);
            }
        }
    }

    // Sort by column if requested
    if (column_sorted) {
        std::ranges::sort(entries, [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.col < b.col) || (a.col == b.col && a.row < b.row);
        });
    }

    return entries;
}

std::vector<MatrixEntry> CSRMatrix::getStrictUpperTriangularPart(const bool column_sorted) const {
    std::vector<MatrixEntry> entries;

    for (int row = 0; row < rows; ++row) {
        for (int idx = this->rowPtr[row]; idx < this->rowPtr[row + 1]; ++idx) {
            if (int col = this->colIndex[idx]; col > row) { // Strict upper triangular condition
                entries.emplace_back(row, col, this->val[idx]);
            }
        }
    }

    // Sort by column if requested
    if (column_sorted) {
        std::ranges::sort(entries, [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.col < b.col) || (a.col == b.col && a.row < b.row);
        });
    }

    return entries;
}

DenseVector CSRMatrix::operator*(const DenseVector &x) const {
    // Check if dimensions match for multiplication
    if (x.size != cols) {
        throw std::invalid_argument("Vector size does not match matrix dimensions.");
    }

    // Result vector initialized with zeros
    DenseVector result(rows);

    // Multiply each row of the matrix by the vector
    for (int row = 0; row < rows; ++row) {
        for (int idx = this->rowPtr[row]; idx < this->rowPtr[row + 1]; ++idx) {
            int col = this->colIndex[idx];
            result.data[row] += this->val[idx] * x.data[col];
        }
    }

    return result;
}

std::vector<MatrixEntry> CSRMatrix::getRow(int row) const{
    if (row < 0 || row >= this->rows) {
        throw std::out_of_range("Row index out of range");
    }

    std::vector<MatrixEntry> rowEntries;
    int start = this->rowPtr[row];
    int end = this->rowPtr[row + 1];

    for (int idx = start; idx < end; ++idx) {
        rowEntries.emplace_back(row, this->colIndex[idx], this->val[idx]);
    }

    return rowEntries;
}


/*
//a_i representing i th row of a and a_k represents k th row of A
//kp1_index is the first index in a_k where the column value of a_k is greater than k
//a_i and a_i_lvl have to have the same size
//a_k and a_k_lvl have to have the same size
//a_k should only contain row values of k
//a_i should only contain row values of i
//a_k is part of the new Matrix
//a_i is under construction
//A_k_lvl is part of the new Matrix
//a_i_lvl is under construction
//new_diagonals is under construction. Has to be initialized with full size already
void CSRMatrix::jloop(int allowed_level,
                      int kp1,
                      std::vector<MatrixEntry> &a_i,
                      std::vector<MatrixEntry> &a_k,
                      std::vector<MatrixLevel> &a_i_lvl,
                      std::vector<MatrixLevel> &a_k_lvl,
                      std::vector<double> &new_diagonals,
                      int a_ik_lvl,
                      double a_ik)
    {
    int current_index = 0;
    const auto end_row = a_k.size();
    int i = a_i[0].row;
    for (int j_idx = 0; j_idx<end_row; ++j_idx) { //todo kp1 vorher schon rausfiltern
        if(a_k[j_idx].col>=kp1) {
            auto& [should_be_k, j, a_kj_val] = a_k[j_idx];

            while(a_i[current_index].col < j && current_index < a_i.size()) {
                ++current_index ;
            }

            if(a_i[current_index].col == j) {
                int LvlEval = std::ranges::min(a_i_lvl[current_index].lvl,1 + a_ik_lvl + a_k_lvl[j_idx].lvl );
                if(LvlEval <= allowed_level) {
                    a_i_lvl[current_index].lvl = LvlEval;
                    a_i[current_index].value -= a_ik*a_kj_val;
                    if (i == j) {
                        new_diagonals[i]=a_i[current_index].value;
                    }
                }
            } else {
                int LvlEval = 1 + a_ik_lvl + a_k_lvl[j_idx].lvl;

                if(LvlEval <= allowed_level) {
                    MatrixEntry tmp(i, j, -a_ik*a_kj_val);
                    auto it = a_i.begin()+ current_index;
                    a_i.insert( it, tmp);
                    MatrixLevel tmplvl(i, j, LvlEval);
                    auto itlvl = a_i_lvl.begin()+ current_index;
                    a_i_lvl.insert(itlvl, tmplvl);
                    if (i == j) {
                        new_diagonals[i]=tmp.value;
                    }
                }

            }
        }


    }
}

//a_i is i th row of A
//new diagonals is initialized with correct size and correct value for first diagonal entry
//do not filter a_i in i loop beause creating new Matrix would be broken;
void CSRMatrix::kloop(int allowed_level, int im1, std::vector<MatrixEntry> a_i,std::vector<MatrixLevel> a_i_lvl,
    std::vector<double> & new_diagonals, std::vector<MatrixEntry> & new_a, std::vector<int>& new_a_row_index,
    std::vector<MatrixLevel>& new_a_lvl) {

    for(int k_idx = 0;  a_i[k_idx].col < im1 && k_idx < a_i.size(); ++k_idx) {

        auto& [ignored, k, a_ik_value] = a_i[k_idx]; //todo check if reference works as intended
        auto& [ignored2, ignored3, a_ik_lvl] = a_i_lvl[k_idx];

        a_i[k_idx].value = a_i[k_idx].value/new_diagonals[k];


        int k_row_index_start = new_a_row_index[k];
        int k_row_index_end = new_a_row_index[k+1];
        std::vector<MatrixEntry> a_k(new_a.begin() + k_row_index_start, new_a.begin() + k_row_index_end);
        std::vector<MatrixLevel> a_k_lvl(new_a_lvl.begin() + k_row_index_start, new_a_lvl.begin() + k_row_index_end);

        jloop(allowed_level,k+1, a_i,
                a_k, a_i_lvl, a_k_lvl, new_diagonals,a_ik_lvl, a_i[k_idx].value );

    }
    new_a.insert(new_a.end(),a_i.begin(),a_i.end());
    new_a_lvl.insert(new_a_lvl.end(),a_i_lvl.begin(),a_i_lvl.end());
    new_a_row_index[im1+1] = new_a.size();
}

std::vector<MatrixEntry> CSRMatrix::iluk(int level) {
//Realisierung des General static ilu patterns

//0. For each (i, j ) ∈ P set aij = 0
//1.For k = 1, . . . , n − 1, Do
//2.  For i = k + 1, n and if (i, k) /∈ P , Do
//3.  aik := aik /a kk
//4.      For j = k + 1, . . . , n and for (i, j ) /∈ P , Do
//5.      a ij := aij − a ik ∗ a kj
//6.      EndDo
//7.  EndDo
//8. EndDo


    std::vector<double> new_diagonals = this->standarddiagonal(); //1. Diagonaleintrag vorbereiten

    std::vector<MatrixEntry> NewMatrix = this->getRow(0); //Neue Matrix mit Reihe 0 Initialisieren
    std::vector<int> NewMatrixRowIndex(this->rows + 1);
    NewMatrixRowIndex[0] = 0;
    NewMatrixRowIndex[1] = NewMatrix.size();
    std::vector<MatrixLevel> matrixLevel;
    for (auto [row, col, val] : NewMatrix) {
        matrixLevel.emplace_back(row, col, 0);
    }

    for (int i = 1; i < rows; ++i) {
        std::vector<MatrixEntry> a_i = this->getRow(i);

        std::vector<MatrixLevel> a_i_lvl;               //lvl für row beziehen
        for (auto [row, col, val] : a_i) {
            a_i_lvl.emplace_back(row, col, 0);
        }


        kloop(level, i, a_i, a_i_lvl, new_diagonals, NewMatrix, NewMatrixRowIndex, matrixLevel);

    }
    return NewMatrix;
}
*/


//a_i representing i th row of a and a_k represents k th row of A
//kp1_index is the first index in a_k where the column value of a_k is greater than k
//a_i and a_i_lvl have to have the same size
//a_k and a_k_lvl have to have the same size
//a_k should only contain row values of k
//a_i should only contain row values of i
//a_k is part of the new Matrix
//a_i is under construction
//A_k_lvl is part of the new Matrix
//a_i_lvl is under construction
//new_diagonals is under construction. Has to be initialized with full size already
void CSRMatrix::jloop(int allowed_level,
                      int kp1,
                      std::vector<MatrixEntry> &a_i,
                      std::vector<MatrixEntry> &a_k,
                      std::vector<MatrixLevel> &a_i_lvl,
                      std::vector<MatrixLevel> &a_k_lvl,
                      std::vector<double> &new_diagonals,
                      int a_ik_lvl,
                      double a_ik)
    {
    int current_index = 0;
    const auto end_row = a_k.size();
    const int i = a_i[0].row;
    for (int j_idx = 0; j_idx<end_row; ++j_idx) { //todo kp1 vorher schon rausfiltern
        if(a_k[j_idx].col>=kp1) {
            auto& [should_be_k, j, a_kj_val] = a_k[j_idx];

            while(current_index < a_i.size() && a_i[current_index].col < j ) {
                ++current_index ;
            }

            if(a_i[current_index].col == j) {
                int LvlEval = std::ranges::min(a_i_lvl[current_index].lvl,1 + a_ik_lvl + a_k_lvl[j_idx].lvl );
                a_i_lvl[current_index].lvl = LvlEval;
                a_i[current_index].value -= a_ik*a_kj_val;
                if (i == j) {
                    new_diagonals[i]=a_i[current_index].value;
                }
            } else {
                int LvlEval = 1 + a_ik_lvl + a_k_lvl[j_idx].lvl;

                MatrixEntry tmp(i, j, -a_ik*a_kj_val);
                a_i.insert( a_i.begin()+ current_index, tmp);
                MatrixLevel tmplvl(i, j, LvlEval);
                a_i_lvl.insert(a_i_lvl.begin()+ current_index, tmplvl);
                if (i == j) {
                    new_diagonals[i]=tmp.value;
                }
            }
        }
    }
}

//a_i is i th row of A
//new diagonals is initialized with correct size and correct value for first diagonal entry
//do not filter a_i in i loop beause creating new Matrix would be broken;
void CSRMatrix::kloop(int allowed_level, int i, std::vector<MatrixEntry> a_i,std::vector<MatrixLevel> a_i_lvl,
    std::vector<double> & new_diagonals, std::vector<MatrixEntry> & new_a, std::vector<int>& new_a_row_index,
    std::vector<MatrixLevel>& new_a_lvl) {

    for(int k_idx = 0;  a_i[k_idx].col < i && k_idx < a_i.size(); ++k_idx) {

        auto& [ignored, k, a_ik_value] = a_i[k_idx];
        auto& [ignored2, ignored3, a_ik_lvl] = a_i_lvl[k_idx];

        if(a_ik_lvl<=allowed_level) {

            a_ik_value = a_ik_value/new_diagonals[k];

            int k_row_index_start = new_a_row_index[k];
            int k_row_index_end = new_a_row_index[k+1];
            std::vector<MatrixEntry> a_k(new_a.begin() + k_row_index_start, new_a.begin() + k_row_index_end);
            std::vector<MatrixLevel> a_k_lvl(new_a_lvl.begin() + k_row_index_start, new_a_lvl.begin() + k_row_index_end);


            jloop(allowed_level,k+1, a_i,
                    a_k, a_i_lvl, a_k_lvl, new_diagonals,a_ik_lvl, a_ik_value );
        }
    }

    //auto startTime = std::chrono::high_resolution_clock::now();
    std::vector<MatrixEntry> filtered_a_i;
    filtered_a_i.reserve(a_i.size());
    std::vector<MatrixLevel> filtered_a_i_lvl;
    filtered_a_i_lvl.reserve(a_i_lvl.size());
    for (int idx = 0; idx< a_i.size(); ++idx) {
        if(a_i_lvl[idx].lvl <= allowed_level) {
            filtered_a_i.push_back(a_i[idx]);
            filtered_a_i_lvl.push_back(a_i_lvl[idx]);
        }
    }
    //auto endTime = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
    //if (elapsed.count() > 0.5) {
    //    std::cout << "slow" << std::endl;
    //}

    new_a.insert(new_a.end(),filtered_a_i.begin(),filtered_a_i.end());
    new_a_lvl.insert(new_a_lvl.end(),filtered_a_i_lvl.begin(),filtered_a_i_lvl.end());
    new_a_row_index[i+1] = new_a.size();
}

std::vector<MatrixEntry> CSRMatrix::iluk(int level) {
//Realisierung des General static ilu patterns

//0. For each (i, j ) ∈ P set aij = 0
//1.For k = 1, . . . , n − 1, Do
//2.  For i = k + 1, n and if (i, k) /∈ P , Do
//3.  aik := aik /a kk
//4.      For j = k + 1, . . . , n and for (i, j ) /∈ P , Do
//5.      a ij := aij − a ik ∗ a kj
//6.      EndDo
//7.  EndDo
//8. EndDo


    std::vector<double> new_diagonals = this->standarddiagonal(); //1. Diagonaleintrag vorbereiten

    std::vector<MatrixEntry> NewMatrix = this->getRow(0); //Neue Matrix mit Reihe 0 Initialisieren
    std::vector<int> NewMatrixRowIndex(this->rows + 1);
    NewMatrixRowIndex[0] = 0;
    NewMatrixRowIndex[1] = NewMatrix.size();
    std::vector<MatrixLevel> matrixLevel;
    for (auto [row, col, val] : NewMatrix) {
        matrixLevel.emplace_back(row, col, 0);
    }

    for (int i = 1; i < rows; ++i) {
        std::vector<MatrixEntry> a_i = this->getRow(i);

        std::vector<MatrixLevel> a_i_lvl;               //lvl für row beziehen
        for (auto [row, col, val] : a_i) {
            a_i_lvl.emplace_back(row, col, 0);
        }

        kloop(level, i, a_i, a_i_lvl, new_diagonals, NewMatrix, NewMatrixRowIndex, matrixLevel);
    }
    return NewMatrix;
}
