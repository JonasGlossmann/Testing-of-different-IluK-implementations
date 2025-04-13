//
// Created by jonas on 27.11.2024.
//

#include "DenseVector.h"

DenseVector::DenseVector()= default;

DenseVector::DenseVector(const int size) {
    this->size = size;
    this->data.resize(this->size,0.0);
}

DenseVector::DenseVector(const int size, const double initValue) {
    this->size = size;
    this->data.resize(this->size,initValue);
}

DenseVector::DenseVector(std::vector<double>input) {
    this->size = input.size();
    this->data = input;
}

DenseVector::DenseVector(const std::string& path) {
    std::ifstream infile(path);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + path);
    }

    std::string line;

    // Skip comments (lines starting with '%')
    while (std::getline(infile, line)) {
        if (line[0] != '%') break;
    }

    // Parse matrix dimensions (for vector, rows = n, columns = 1)
    std::istringstream iss(line);
    int rows, cols, nnz;  // nnz is the number of non-zero elements, for a vector this is rows (assuming dense)
    iss >> rows >> cols >> nnz;

    // Initialize a vector to hold the values
    this->data.resize(rows, 0.0);  // Dense vector initialized with zeros

    for (int i = 0; i < rows; ++i) {
        if (!(infile >> this->data[i])) {
            throw std::runtime_error("Error reading vector values from file.");
        }
    }

    this->size = rows;
}

double DenseVector::norm() {
    double sum = 0.0;
    for (const double val : this->data) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

double DenseVector::dot(const DenseVector& other_vector) {
    double result = 0.0;
    for (int i = 0; i < this->size; i++) {
        result += this->data[i] * other_vector.data[i];
    }
    return result;
}

DenseVector DenseVector::operator+(const DenseVector &other_vector) {
    DenseVector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.data[i] = this->data[i] + other_vector.data[i];
    }
    return result;
}

DenseVector DenseVector::operator-(const DenseVector &other_vector) {
    DenseVector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.data[i] = this->data[i] - other_vector.data[i];
    }
    return result;
}

DenseVector DenseVector::operator*(const DenseVector &other_vector) {
    DenseVector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.data[i] = this->data[i] * other_vector.data[i];
    }
    return result;
}

/*
void writeDenseVectorToMTX(const Eigen::VectorXd& vector, const std::string& filePath) {
    std::ofstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return;
    }

    // Write the Matrix Market header
    file << "%%MatrixMarket matrix coordinate real general\n";
    file << vector.size() << " 1 " << vector.size() << "\n";

    // Write the vector elements as (row, col, value)
    for (int i = 0; i < vector.size(); ++i) {
        file << (i + 1) << " 1 " << std::setprecision(10) << vector(i) << "\n";
    }

    file.close();
}
*/


