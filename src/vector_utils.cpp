#include "elda/vector_utils.hpp"
#include <vector>

namespace linalg {
    bool check_similar(matrix mat1, matrix mat2) {
        // Ignore zero rows before comparing the remaining row sequence.
        std::vector<double> zero (mat1.col,0);
        std::vector<std::vector<double>> check1;
        for (int i = 0; i < mat1.row; i++) {
            if (mat1.arr[i] != zero) {
                check1.push_back(mat1.arr[i]);
            }
        }
        std::vector<std::vector<double>> check2;
        for (int i = 0; i < mat2.row; i++) {
            if (mat2.arr[i] != zero) {
                check2.push_back(mat2.arr[i]);
            }
        }
        if (check1==check2) {
            return true;
        }
        return false;
    }
}
