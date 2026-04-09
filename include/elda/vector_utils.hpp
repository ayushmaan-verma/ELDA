#ifndef ELDA_VECTOR_UTILS_H
#define ELDA_VECTOR_UTILS_H

#include "matrix.hpp"

namespace linalg {
    /// Compares matrices after removing all-zero rows from each one.
    /// This is a library-specific helper, not the standard notion of matrix similarity.
    bool check_similar(matrix mat1, matrix mat2);
}

#endif // ELDA_VECTOR_UTILS_H
