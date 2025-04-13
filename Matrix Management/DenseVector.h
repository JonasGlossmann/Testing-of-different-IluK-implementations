//
// Created by jonas on 27.11.2024.
//

#ifndef DENSEVECTOR_H
#define DENSEVECTOR_H
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>


class DenseVector {
public:
    std::vector<double> data;
    int size{};

    DenseVector();
    explicit DenseVector(int size);
    DenseVector(int size, double initValue);
    explicit DenseVector(std::vector<double> input);
    explicit DenseVector(const std::string& path);

    double norm();
    DenseVector operator+(const DenseVector& other_vector);
    DenseVector operator-(const DenseVector& other_vector);
    DenseVector friend operator*(double scalar, DenseVector vector) {
        DenseVector result(vector.size);
        for (int i = 0; i < vector.size; i++) {
            result.data[i] = vector.data[i] * scalar;
        }
        return result;
    };
    DenseVector friend operator/(double scalar, DenseVector vector) {
        DenseVector result(vector.size);
        for (int i = 0; i < vector.size; i++) {
            result.data[i] = scalar /vector.data[i] ;
        }
        return result;
    };
    DenseVector operator*(const DenseVector& other_vector);
    double dot(const DenseVector& other_vector);



};



#endif //DENSEVECTOR_H
