//
// Created by jonas on 08.12.2024.
//

#include "CSCMatrix.h"

void CSCMatrix::fromMatrixEntries(const std::vector<MatrixEntry>& matrix, const int cols) {
    
    this->colPtr.resize(cols + 1, 0);  // Start with zero column counts
    this->rowIndex.reserve(matrix.size());
    this->val.reserve(matrix.size());

    // Populate CSC structure
    int current_col = -1; // Tracks the current column being processed
    for (const auto& [row, col, value] : matrix) {
        // Fill in colPtr when moving to a new column
        while (current_col < col) {
            current_col++;
            this->colPtr[current_col] = this->rowIndex.size();
        }

        // Add the row index and value for this entry
        this->rowIndex.push_back(row);
        this->val.push_back(value);
    }

    // Complete colPtr for remaining columns
    while (current_col < cols) {
        current_col++;
        colPtr[current_col] = this->rowIndex.size();
    }
}

CSCMatrix::CSCMatrix(const std::string &filepath) {
    const COOMatrixSorted matrix(filepath);
    this->rows = matrix.rows;
    this->cols = matrix.cols;
    fromMatrixEntries(matrix.getFullMatrix(true), this->cols);
}

CSCMatrix::CSCMatrix(const COOMatrixSorted& matrix) {
    this->rows = matrix.rows;
    this->cols = matrix.cols;
    fromMatrixEntries(matrix.getFullMatrix(true), this->cols);
}

CSCMatrix::CSCMatrix(const CSRMatrix& matrix) {
    const COOMatrixSorted new_matrix(matrix);
    this->rows = new_matrix.rows;
    this->cols = new_matrix.cols;
    fromMatrixEntries(new_matrix.getFullMatrix(true), this->cols);
}

CSCMatrix::CSCMatrix(std::vector<MatrixEntry> matrix, const int rows, const int cols) {
    this->rows = rows;
    this->cols = cols;
    std::ranges::sort(matrix, [](const MatrixEntry& a, const MatrixEntry& b) {
    return (a.col < b.col) || (a.col == b.col && a.row < b.row);
    });
    fromMatrixEntries(matrix, this->cols);
}

DenseVector CSCMatrix::diagonal() {
    DenseVector result(this->cols);
    for (int col = 0; col < this->cols; ++col) {
        for (int idx = this->colPtr[col]; idx < this->colPtr[col + 1]; ++idx) {
            if (this->rowIndex[idx] == col) { // Check for diagonal element
                result.data[col] = this->val[idx];
                break; // No need to check further in this column
            }
        }
    }
    return result;
}

double CSCMatrix::coeffRef(const int row, const int col) {
    // Check column bounds
    if (col < 0 || col >= this->colPtr.size() - 1) {
        throw std::out_of_range("Column index out of range");
    }

    // Get the start and end of the column in the CSC structure
    int start = this->colPtr[col];
    int end = this->colPtr[col + 1];

    // Search for the row in this column
    for (int idx = start; idx < end; ++idx) {
        if (this->rowIndex[idx] == row) {
            return this->val[idx]; // Return the value if the row matches
        }
    }

    // Return 0.0 if no match is found (sparse assumption)
    return 0.0;
}

std::vector<MatrixEntry> CSCMatrix::getFullMatrix(const bool column_sorted) const {
    std::vector<MatrixEntry> result;

    for (int col = 0; col < cols; ++col) {
        // Iterate through the non-zero elements in this column
        for (int idx = this->colPtr[col]; idx < this->colPtr[col + 1]; ++idx) {
            result.emplace_back(this->rowIndex[idx], col, this->val[idx]);
        }
    }

    if (!column_sorted) {
        std::ranges::sort(result, [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.row < b.row) || (a.row == b.row && a.col < b.col);
        });
    }

    return result;
}

std::vector<MatrixEntry> CSCMatrix::getStrictUpperTriangularPart(const bool column_sorted) const {
    std::vector<MatrixEntry> result = this->getFullMatrix(column_sorted);

    auto it = std::ranges::remove_if(result, [](const MatrixEntry& entry) {
        return entry.row >= entry.col; // Remove if not in strict upper triangular part
    }).begin();

    // Erase the "removed" entries from the vector
    result.erase(it, result.end());
    if (!column_sorted) {
        std::ranges::sort(result, [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.row < b.row) || (a.row == b.row && a.col < b.col);
        });
    }
    return result;
}

std::vector<MatrixEntry> CSCMatrix::getStrictLowerTriangularPart(const bool column_sorted) const {
    std::vector<MatrixEntry> result = this->getFullMatrix(column_sorted);

    const auto it = std::ranges::remove_if(result, [](const MatrixEntry& entry) {
        return entry.row <= entry.col; // Remove if not in strict upper triangular part
    }).begin();

    // Erase the "removed" entries from the vector
    result.erase(it, result.end());
    if (!column_sorted) {
        std::ranges::sort(result, [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.row < b.row) || (a.row == b.row && a.col < b.col);
        });
    }
    return result;
}

DenseVector CSCMatrix::operator*(const DenseVector &x) const {
    DenseVector result(this->cols);

    for (int col = 0; col < this->cols; ++col) {
        for (int idx = this->colPtr[col]; idx < this->colPtr[col + 1]; ++idx) {
            result.data[this->rowIndex[idx]] += this->val[idx] * x.data[col];
        }
    }
    return result;
}
