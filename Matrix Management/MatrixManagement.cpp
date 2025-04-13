//
// Created by jonas on 07.11.2024.
//

#include "MatrixManagement.h"

namespace fs=std::filesystem;

std::string wstring_to_string(const std::wstring& wstr) {
    // Use wstring_convert to convert wstring to string
    std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
    return converter.to_bytes(wstr);
}


MatrixManagement::MatrixManagement(const std::wstring& root_dir,
                     const std::wstring& matrix_name,
                     const std::wstring& rhs_name,
                     const std::wstring& startvector_name) {
    this->currentRow = 0;

    // Iterate over directories and subdirectories
    for (const auto& dir_entry : fs::recursive_directory_iterator(root_dir)) {
        if (fs::is_directory(dir_entry)) {
            bool f1_found = false, f2_found = false, f3_found = false;
            std::wstring f1_path, f2_path, f3_path;

            // Check each file within the current directory
            for (const auto& file_entry : fs::directory_iterator(dir_entry)) {
                if (file_entry.path().filename() == matrix_name) {
                    f1_found = true;
                    f1_path = file_entry.path().wstring();
                } else if (file_entry.path().filename() == rhs_name) {
                    f2_found = true;
                    f2_path = file_entry.path().wstring();
                } else if (file_entry.path().filename() == startvector_name) {
                    f3_found = true;
                    f3_path = file_entry.path().wstring();
                }
            }

            // If all three files are found, add them as a triple
            if (f1_found && f2_found && f3_found) {
                matrices.push_back({f1_path, f2_path, f3_path});
            }
        }
    }

}

bool MatrixManagement::validIndex() const {
    return this->matrices.size() >= this->currentRow ;
}

void MatrixManagement::next() {
    currentRow++;
}

MatrixTriple MatrixManagement::getMatrix() {
    auto [matrix, rhs, startvector_name] = this->matrices[this->currentRow];
    MatrixTriple triple(wstring_to_string(matrix));

    triple.rhs = DenseVector(wstring_to_string(rhs));
    triple.startvector = DenseVector(wstring_to_string(startvector_name));
    return triple;
}









