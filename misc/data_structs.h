//
// Created by jonas on 05.04.2025.
//

#ifndef DATA_STRUCTS_H
#define DATA_STRUCTS_H

struct MatrixEntry {
    int row;
    int col;
    double value;
};

struct MatrixLevel {
    int row;
    int col;
    int lvl;
};

struct FileTriple {
    std::wstring matrix ;
    std::wstring rhs;
    std::wstring startvector;
};


#endif //DATA_STRUCTS_H
