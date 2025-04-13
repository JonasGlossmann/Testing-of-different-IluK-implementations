//
// Created by jonas on 27.11.2024.
//

#include "COOMatrixSorted.h"

#include "../misc/TimeControl.h"


struct PairHash {
    std::size_t operator()(const std::pair<int, int>& p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

double COOMatrixSorted::binary_search(const int target_row, const int target_col) {
    MatrixEntry target = {target_row, target_col, 0.0};  // The value is irrelevant for comparison

    // Perform binary search
    auto it = std::lower_bound(this->matrix.begin(), this->matrix.end(), target,
        [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.col < b.col) || (a.col == b.col && a.row < b.row);
        });

    // Check if the found entry matches the target row and column
    if (it != this->matrix.end() && it->row == target_row && it->col == target_col) {
        return it->value;
    }

    // Return empty optional if not found
    return 0.0;
}

void COOMatrixSorted::get_diagonal_pointer_pairs() {

    for (auto&[row, col, val] : this->matrix) {
        if (row == col) {
            this->diags.emplace_back(row, &val);
        }
    }
}

//read the matrix from a .mtx file
void COOMatrixSorted::readFromMTX(const std::string& filepath) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filepath);
    }

    std::string line;
    // Skip comments
    while (std::getline(infile, line)) {
        if (line[0] != '%') break;
    }

    // Read matrix dimensions and number of non-zero elements
    std::istringstream iss(line);
    int nnz; // number of non-zero elements
    iss >> this->rows >> this->cols >> nnz;

    // Read the non-zero values and their positions
    int row, col;
    double value;
    std::vector<MatrixEntry> temp;
    this->matrix.reserve(nnz);
    while (infile >> row >> col >> value) {
        this->matrix.push_back({row-1, col-1, value});
    }
    std::ranges::sort(this->matrix, [](const MatrixEntry& a, const MatrixEntry& b) {
    return (a.col < b.col) || (a.col == b.col && a.row < b.row);
    });
}

COOMatrixSorted::COOMatrixSorted(const std::string &filepath) {
    this->readFromMTX(filepath);
    this->get_diagonal_pointer_pairs();
}

COOMatrixSorted::COOMatrixSorted(const CSRMatrix& matrix) {
    this->rows = matrix.rows;
    this->cols = matrix.cols;
    this->matrix = matrix.getFullMatrix(true);
    this->get_diagonal_pointer_pairs();
}

COOMatrixSorted::COOMatrixSorted(const CSCMatrix& matrix) {
    this->rows = matrix.rows;
    this->cols = matrix.cols;
    this->matrix = matrix.getFullMatrix(true);
    this->get_diagonal_pointer_pairs();
}

COOMatrixSorted::COOMatrixSorted(const std::vector<MatrixEntry>& matrix, const int rows, const int cols,
                    const bool column_sorted) {
    this->rows = rows;
    this->cols = cols;
    this->matrix = matrix;
    if(column_sorted) {
        std::ranges::sort(this->matrix, [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.col < b.col) || (a.col == b.col && a.row < b.row);
        });
    } else {
        std::ranges::sort(this->matrix, [](const MatrixEntry& a, const MatrixEntry& b) {
        return (a.row < b.row) || (a.row == b.row && a.col < b.col);
        });
    }
    this->get_diagonal_pointer_pairs();
}

// Matrix-vector multiplication: y = A * x
DenseVector COOMatrixSorted::operator*(const DenseVector &x) const {
    if (x.size != this->cols) {
        std::cout<< x.size << "    " <<cols<< std::endl;
        throw std::invalid_argument("Vector size must match number of columns in the matrix.");
    }

    DenseVector y(rows); // Result vector initialized to zero

    for (auto [row, col,val] : this->matrix) {
        y.data[row] += x.data[col] * val;
    }

    return y;
}

DenseVector COOMatrixSorted::diagonal() {
    auto x = DenseVector(this->rows);
    for (auto [idx, val] : this->diags) {
        x.data[idx] = val[0];
    }
    return x;
}

double COOMatrixSorted::coeffRef(const int row, const int col) {
    return this->binary_search(row, col);
}

std::vector<MatrixEntry> COOMatrixSorted::getFullMatrix(const bool column_sorted) const {

    if (column_sorted) {
        return this->matrix;
    }
    // sorting indexing for row format
    std::vector<MatrixEntry> Matrix = this->matrix;

    std::ranges::sort(Matrix, [](const MatrixEntry& a, const MatrixEntry& b) {
    return (a.row < b.row) || (a.row == b.row && a.col < b.col);
    });
    return Matrix;
}

std::vector<MatrixEntry> COOMatrixSorted::getStrictUpperTriangularPart(const bool column_sorted) const {

    std::vector<MatrixEntry> Matrix;
    for (auto [row, col, val] : this->matrix) {
        if(row < col) {
            Matrix.push_back({row,col,val});
        }
    }

    if (column_sorted) {
        return Matrix;
    }

    // sorting indexing for row format

    std::ranges::sort(Matrix, [](const MatrixEntry& a, const MatrixEntry& b) {
    return (a.row < b.row) || (a.row == b.row && a.col < b.col);
    });
    return Matrix;
}

std::vector<MatrixEntry> COOMatrixSorted::getStrictLowerTriangularPart(const bool column_sorted) const {
    std::vector<MatrixEntry> Matrix;
    for (auto [row, col, val] : this->matrix) {
        if(row > col) {
            Matrix.push_back({row,col,val});
        }
    }

    if (column_sorted) {
        return Matrix;
    }

    // sorting indexing for row format

    std::ranges::sort(Matrix, [](const MatrixEntry& a, const MatrixEntry& b) {
    return (a.row < b.row) || (a.row == b.row && a.col < b.col);
    });
    return Matrix;
}

