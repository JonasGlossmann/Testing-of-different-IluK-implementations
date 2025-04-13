//
// Created by jonas on 07.11.2024.
//

#ifndef MATRIXMANAGEMENT_H
#define MATRIXMANAGEMENT_H
#include <filesystem>
#include <string>
#include <vector>

#include "DenseVector.h"
#include "../misc/data_structs.h"
#include "COOMatrixSorted.h"



class MatrixManagement {
    std::vector<FileTriple> matrices;
    int currentRow = 0;

    public:
    MatrixManagement(const std::wstring& root_dir,
                     const std::wstring& matrix_name,
                     const std::wstring& rhs_name,
                     const std::wstring& startvector_name);
    bool validIndex() const;
    void next();
    MatrixTriple getMatrix();

};


#endif //MATRIXMANAGEMENT_H
