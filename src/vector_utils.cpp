#include "elda/vector_utils.hpp"
#include <stdexcept>

namespace linalg {

// Convenience constructors for zero-initialized column vectors.
matrix vec1() { return matrix(1, 1); }
matrix vec2() { return matrix(2, 1); }
matrix vec3() { return matrix(3, 1); }
matrix vec4() { return matrix(4, 1); }
matrix vec5() { return matrix(5, 1); }

// Convenience constructors for value-initialized column vectors.
matrix vec1(double x) {
    matrix m(1, 1);
    m.arr[0] = x;
    return m;
}

matrix vec2(double x, double y) {
    matrix m(2, 1);
    m.arr[0] = x;
    m.arr[1] = y;
    return m;
}

matrix vec3(double x, double y, double z) {
    matrix m(3, 1);
    m.arr[0] = x;
    m.arr[1] = y;
    m.arr[2] = z;
    return m;
}

matrix vec4(double x, double y, double z, double w) {
    matrix m(4, 1);
    m.arr[0] = x;
    m.arr[1] = y;
    m.arr[2] = z;
    m.arr[3] = w;
    return m;
}

matrix vec5(double x, double y, double z, double w, double u) {
    matrix m(5, 1);
    m.arr[0] = x;
    m.arr[1] = y;
    m.arr[2] = z;
    m.arr[3] = w;
    m.arr[4] = u;
    return m;
}

bool check_lin_comb(matrix m1, matrix m2) {
    if (m1.get_cols() != 1 || m2.get_cols() != 1) {
        throw std::runtime_error("All inputs must be column vectors (exactly 1 column).");
    }
    if (m1.get_rows() != m2.get_rows()) {
        throw std::runtime_error("Dimension mismatch: All vectors must have the same number of rows.");
    }
    matrix aug(m1.get_rows(), 2);
    for (size_t i = 0; i < m1.get_rows(); i++) {
        aug.arr[i * aug.col] = m1.arr[i];
        aug.arr[i * aug.col + 1] = m2.arr[i];
    }
    matrix coeff(m1.get_rows(), 1);
    for (size_t i = 0; i < m1.get_rows(); i++) {
        coeff.arr[i * coeff.col] = m1.arr[i];
    }
    return aug.rank() == coeff.rank();
}

bool check_lin_comb(matrix m1, matrix m2, matrix m3) {
    if (m1.get_cols() != 1 || m2.get_cols() != 1 || m3.get_cols() != 1) {
        throw std::runtime_error("All inputs must be column vectors (exactly 1 column).");
    }
    if (m1.get_rows() != m2.get_rows() || m2.get_rows() != m3.get_rows()) {
        throw std::runtime_error("Dimension mismatch: All vectors must have the same number of rows.");
    }
    matrix aug(m1.get_rows(), 3);
    for (size_t i = 0; i < m1.get_rows(); i++) {
        aug.arr[i * aug.col] = m1.arr[i];
        aug.arr[i * aug.col + 1] = m2.arr[i];
        aug.arr[i * aug.col + 2] = m3.arr[i];
    }
    matrix coeff(m1.get_rows(), 2);
    for (size_t i = 0; i < m1.get_rows(); i++) {
        coeff.arr[i * coeff.col] = m1.arr[i];
        coeff.arr[i * coeff.col + 1] = m2.arr[i];
    }
    return aug.rank() == coeff.rank();
}

bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4) {
    if (m1.get_cols() != 1 || m2.get_cols() != 1 || m3.get_cols() != 1 || m4.get_cols() != 1) {
        throw std::runtime_error("All inputs must be column vectors (exactly 1 column).");
    }
    if (m1.get_rows() != m2.get_rows() || m2.get_rows() != m3.get_rows() || m3.get_rows() != m4.get_rows()) {
        throw std::runtime_error("Dimension mismatch: All vectors must have the same number of rows.");
    }
    matrix aug(m1.get_rows(), 4);
    for (size_t i = 0; i < m1.get_rows(); i++) {
        aug.arr[i * aug.col] = m1.arr[i];
        aug.arr[i * aug.col + 1] = m2.arr[i];
        aug.arr[i * aug.col + 2] = m3.arr[i];
        aug.arr[i * aug.col + 3] = m4.arr[i];
    }
    matrix coeff(m1.get_rows(), 3);
    for (size_t i = 0; i < m1.get_rows(); i++) {
        coeff.arr[i * coeff.col] = m1.arr[i];
        coeff.arr[i * coeff.col + 1] = m2.arr[i];
        coeff.arr[i * coeff.col + 2] = m3.arr[i];
    }
    return aug.rank() == coeff.rank();
}

bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4, matrix m5) {
    if (m1.get_cols() != 1 || m2.get_cols() != 1 || m3.get_cols() != 1 || m4.get_cols() != 1 || m5.get_cols() != 1) {
        throw std::runtime_error("All inputs must be column vectors (exactly 1 column).");
    }
    if (m1.get_rows() != m2.get_rows() || m2.get_rows() != m3.get_rows() || m3.get_rows() != m4.get_rows() || m4.get_rows() != m5.get_rows()) {
        throw std::runtime_error("Dimension mismatch: All vectors must have the same number of rows.");
    }
    matrix aug(m1.get_rows(), 5);
    for (size_t i = 0; i < m1.get_rows(); i++) {
        aug.arr[i * aug.col] = m1.arr[i];
        aug.arr[i * aug.col + 1] = m2.arr[i];
        aug.arr[i * aug.col + 2] = m3.arr[i];
        aug.arr[i * aug.col + 3] = m4.arr[i];
        aug.arr[i * aug.col + 4] = m5.arr[i];
    }
    matrix coeff(m1.get_rows(), 4);
    for (size_t i = 0; i < m1.get_rows(); i++) {
        coeff.arr[i * coeff.col] = m1.arr[i];
        coeff.arr[i * coeff.col + 1] = m2.arr[i];
        coeff.arr[i * coeff.col + 2] = m3.arr[i];
        coeff.arr[i * coeff.col + 3] = m4.arr[i];
    }
    return aug.rank() == coeff.rank();
}

bool check_lin_comb(matrix target, const std::vector<matrix>& vectors) {
    if (target.get_cols() != 1) {
        throw std::runtime_error("Target matrix must be a column vector (exactly 1 column).");
    }
    for (const auto& vec : vectors) {
        if (vec.get_cols() != 1) {
            throw std::runtime_error("All input matrices in the vector set must be column vectors (exactly 1 column).");
        }
        if (vec.get_rows() != target.get_rows()) {
            throw std::runtime_error("Dimension mismatch: All vectors must have the same number of rows.");
        }
    }
    
    if (vectors.empty()) return false;
    
    size_t rows = target.get_rows();
    size_t cols = vectors.size();
    
    matrix coeff(rows, cols);
    matrix aug(rows, cols + 1);
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            coeff.arr[i * coeff.col + j] = vectors[j].arr[i];
            aug.arr[i * aug.col + j] = vectors[j].arr[i];
        }
        aug.arr[i * aug.col + cols] = target.arr[i];
    }
    return aug.rank() == coeff.rank();
}

} // namespace linalg