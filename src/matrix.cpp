#include "elda/matrix.hpp"
#include <algorithm>
#include <string>
#include <sstream>
#include <tuple>
#include <stdexcept>
#include <cmath>
#include <ostream>

namespace linalg {

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------

void matrix::input() {
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            std::cin >> arr[i * col + j];
}

void matrix::print() {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++)
            std::cout << arr[i * col + j] << " ";
        std::cout << std::endl;
    }
}

std::string matrix::to_string() const {
    std::ostringstream oss;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            oss << arr[i * col + j];
            if (j < col - 1) oss << " ";
        }
        oss << "\n";
    }
    return oss.str();
}

std::ostream& operator<<(std::ostream& os, const matrix& m) {
    for (int i = 0; i < m.row; i++) {
        for (int j = 0; j < m.col; j++)
            os << m(static_cast<size_t>(i), static_cast<size_t>(j)) << " ";
        os << std::endl;
    }
    return os;
}

// ---------------------------------------------------------------------------
// Utility free functions
// ---------------------------------------------------------------------------

matrix identity(int n) {
    if (n < 0)
        throw std::runtime_error("Identity matrix size must be non-negative.");
    matrix m(n, n);
    for (int i = 0; i < n; i++)
        m.arr[i * n + i] = 1.0;
    return m;
}

void neg_zero(matrix& m) {
    for (auto& v : m.arr)
        if (v == -0.0) v = 0.0;
}

void fpg(matrix& m) {
    for (auto& v : m.arr)
        if (std::abs(v) <= EPS) v = 0.0;
}

// ---------------------------------------------------------------------------
// Approximate-comparison helpers
// ---------------------------------------------------------------------------

/// Returns true when m1 and m2 have the same shape and every element pair
/// satisfies |m1[i] - m2[i]| <= tol.  Shape mismatches return false without
/// accessing any element (safe for empty or zero-sized matrices).
bool matrices_approx_equal(const matrix& m1, const matrix& m2, double tol) {
    if (m1.row != m2.row || m1.col != m2.col) return false;
    const size_t n = static_cast<size_t>(m1.row) * static_cast<size_t>(m1.col);
    for (size_t i = 0; i < n; ++i)
        if (!approx_equal(m1.arr[i], m2.arr[i], tol)) return false;
    return true;
}

/// Exact equality: delegates to matrices_approx_equal with zero tolerance so
/// shape mismatches are still caught safely before any element access.
matrix matrix::operator=(const matrix& m2) {
    if ((row != m2.row) || (col != m2.col)) {
        throw std::runtime_error("Assignment Operator : Dimension Mismatch in LHS & RHS !!!");
    }
    arr = m2.arr;
    return *this;
}

matrix matrix::operator+(const matrix& m2) {
    if ((row != m2.row) || (col != m2.col)) {
        throw std::runtime_error("Addition defined only if both matrices share identical dimensions.");
    }
    matrix m3(row, col);
    size_t total_elements = get_rows() * get_cols();
    // Cache locality optimization: Linear element loop pass
    for (size_t i = 0; i < total_elements; i++) {
        m3.arr[i] = this->arr[i] + m2.arr[i];
    }
    return m3;
}

matrix matrix::operator-(const matrix& m2) {
    if ((row != m2.row) || (col != m2.col)) {
        throw std::runtime_error("Subtraction defined only if both matrices share identical dimensions.");
    }
    matrix m3(row, col);
    size_t total_elements = get_rows() * get_cols();
    for (size_t i = 0; i < total_elements; i++) {
        m3.arr[i] = this->arr[i] - m2.arr[i];
    }
    return m3;
}

matrix matrix::operator*(const matrix& m2) {
    if (col != m2.row) {
        throw std::runtime_error("Matrix multiplication dimensions incompatible.");
    }
    matrix m3(row, m2.col);
    // Highly optimized cache-friendly index scan loops
    for (int i = 0; i < row; i++) {
        for (int k = 0; k < col; k++) {
            double r = arr[i * col + k];
            for (int j = 0; j < m2.col; j++) {
                m3.arr[i * m2.col + j] += r * m2.arr[k * m2.col + j];
            }
        }
    }
    fpg(m3);
    return m3;
}

matrix matrix::operator*(double d) {
    matrix m3(*this);
    size_t total_elements = get_rows() * get_cols();
    for (size_t i = 0; i < total_elements; i++) {
        m3.arr[i] *= d;
    }
    fpg(m3);
    return m3;
}

matrix matrix::operator/(double scalar) {
    if (std::abs(scalar) <= EPS)
        throw std::runtime_error("Division by zero in matrix scalar division.");
    return (*this) * (1.0 / scalar);
}

bool operator==(const matrix& m1, const matrix& m2) {
    return matrices_approx_equal(m1, m2, 0.0);
}

bool shape_comp(const matrix& m1, const matrix& m2) {
    return (m1.row == m2.row) && (m1.col == m2.col);
}

// ---------------------------------------------------------------------------
// Row/column helpers
// ---------------------------------------------------------------------------

matrix matrix::get_row(int r) {
    matrix rw(row, col);
    for (int j = 0; j < col; j++)
        rw.arr[r * col + j] = arr[r * col + j];
    return rw;
}

matrix matrix::get_col(int c) {
    matrix cn(row, col);
    for (int i = 0; i < row; i++)
        cn.arr[i * col + c] = arr[i * col + c];
    return cn;
}

matrix matrix::get_row_vec(int r) {
    matrix rw(1, col);
    for (int j = 0; j < col; j++)
        rw.arr[j] = arr[r * col + j];
    return rw;
}

matrix matrix::get_col_vec(int c) {
    matrix cn(row, 1);
    for (int i = 0; i < row; i++)
        cn.arr[i] = arr[i * col + c];
    return cn;
}

void matrix::replace_row(int r, const matrix& rw) {
    for (int j = 0; j < col; j++)
        arr[r * col + j] = rw.arr[j];
}

void matrix::replace_col(int c, const matrix& cn) {
    for (int i = 0; i < row; i++)
        arr[i * col + c] = cn.arr[i];
}

matrix matrix::row_op(int i, double coeff, int j) {
    matrix m(*this);
    for (int k = 0; k < col; k++)
        m.arr[i * col + k] += coeff * arr[j * col + k];
    fpg(m);
    return m;
}

matrix matrix::col_op(int i, double coeff, int j) {
    matrix m(*this);
    for (int k = 0; k < row; k++)
        m.arr[k * col + i] += coeff * m.arr[k * col + j];
    fpg(m);
    return m;
}

matrix matrix::row_swap(int i, int j) {
    matrix m(*this);
    for (int k = 0; k < col; k++)
        std::swap(m.arr[i * col + k], m.arr[j * col + k]);
    fpg(m);
    return m;
}

matrix matrix::col_swap(int i, int j) {
    matrix m(*this);
    for (int k = 0; k < row; k++)
        std::swap(m.arr[k * col + i], m.arr[k * col + j]);
    fpg(m);
    return m;
}

matrix matrix::row_multi(int i, double factor) {
    matrix m(*this);
    for (int k = 0; k < col; k++)
        m.arr[i * col + k] *= factor;
    fpg(m);
    return m;
}

matrix matrix::col_multi(int i, double factor) {
    matrix m(*this);
    for (int k = 0; k < row; k++)
        m.arr[k * col + i] *= factor;
    fpg(m);
    return m;
}

// ---------------------------------------------------------------------------
// Elimination algorithms
// ---------------------------------------------------------------------------

int matrix::echelon() {
    int swaps = 1;
    int r = 0;
    for (int c = 0; c < col && r < row; c++) {
        // Partial pivoting: find largest magnitude in this column
        int pivot_row = r;
        double max_val = std::abs(arr[r * col + c]);
        for (int i = r + 1; i < row; i++) {
            if (std::abs(arr[i * col + c]) > max_val) {
                max_val = std::abs(arr[i * col + c]);
                pivot_row = i;
            }
        }
        if (max_val < 1e-9) continue; // Zero column – skip
        if (pivot_row != r) {
            *this = row_swap(r, pivot_row);
            swaps *= -1;
        }
        for (int i = r + 1; i < row; i++) {
            const double multiplier = arr[i * col + c] / arr[r * col + c];
            *this = row_op(i, -multiplier, r);
        }
        r++;
    }
    neg_zero(*this);
    fpg(*this);
    return swaps;
}

int matrix::gaussian() {
    int swaps = 1;
    int r = 0;
    for (int c = 0; c < col && r < row; c++) {
        int pivot_row = r;
        double max_val = std::abs(arr[r * col + c]);
        for (int i = r + 1; i < row; i++) {
            if (std::abs(arr[i * col + c]) > max_val) {
                max_val = std::abs(arr[i * col + c]);
                pivot_row = i;
            }
        }
        if (max_val < 1e-9) continue;
        if (pivot_row != r) {
            *this = row_swap(r, pivot_row);
            swaps *= -1;
        }
        // Normalize pivot row
        const double divisor = arr[r * col + c];
        for (int k = 0; k < col; k++)
            arr[r * col + k] /= divisor;
        // Eliminate below
        for (int i = r + 1; i < row; i++) {
            const double multiplier = arr[i * col + c];
            *this = row_op(i, -multiplier, r);
        }
        r++;
    }
    neg_zero(*this);
    fpg(*this);
    return swaps;
}

int matrix::gauss_jordan() {
    const int swaps = this->gaussian();
    // Backward elimination from the bottom up
    for (int i = row - 1; i > 0; i--) {
        int pivot_col = -1;
        for (int j = 0; j < col; j++) {
            if (std::abs(arr[i * col + j]) > 1e-9) { pivot_col = j; break; }
        }
        if (pivot_col == -1) continue;
        for (int k = i - 1; k >= 0; k--) {
            const double multiplier = arr[k * col + pivot_col];
            *this = row_op(k, -multiplier, i);
        }
    }
    fpg(*this);
    return swaps;
}

int matrix::canonical() {
    const int swaps = this->gaussian();
    int r = 0;
    for (int c = 0; c < col && r < row; c++) {
        if (std::abs(arr[r * col + c]) > 1e-9) {
            for (int j = c + 1; j < col; j++) {
                const double multiplier = arr[r * col + j];
                *this = col_op(j, -multiplier, c);
            }
            r++;
        }
    }
    fpg(*this);
    return swaps;
}

int matrix::rank() {
    int rank_cnt = 0;
    matrix m(*this);
    m.gaussian();
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            if (std::abs(m.arr[i * col + j]) > 1e-9) { rank_cnt++; break; }
    return rank_cnt;
}

int matrix::free_variable() {
    if (col != row)
        throw std::runtime_error("The dimension of Coefficient Matrix must be N x N.");
    return row - this->rank();
}

// ---------------------------------------------------------------------------
// Scalar properties
// ---------------------------------------------------------------------------

double matrix::trace() {
    if (row != col)
        throw std::runtime_error("Trace is defined only for square matrix.");
    double sum = 0;
    for (int i = 0; i < row; i++) sum += arr[i * col + i];
    return sum;
}

double matrix::diag_product() {
    double prod = 1.0;
    for (int i = 0; i < std::min(row, col); i++)
        prod *= arr[i * col + i];
    return prod;
}

double matrix::det() {
    if (row != col)
        throw std::runtime_error("Determinant is defined only for square matrix.");
    matrix m(*this);
    double determinant = m.echelon();
    determinant *= m.diag_product();
    if (determinant == -0.0) determinant = 0.0;
    return determinant;
}

double matrix::minor(int x, int y) {
    if (row != col)
        throw std::runtime_error("Cofactor & Minor defined only for square matrix.");
    matrix mat(row - 1, col - 1);
    int rr = 0;
    for (int i = 0; i < row; i++) {
        if (i == x) continue;
        int cc = 0;
        for (int j = 0; j < col; j++) {
            if (j == y) continue;
            mat.arr[rr * (col - 1) + cc] = arr[i * col + j];
            cc++;
        }
        rr++;
    }
    return mat.det();
}

double matrix::cofactor(int x, int y) {
    int check = x + y + 2;
    return (check % 2) ? -this->minor(x, y) : this->minor(x, y);
}

// ---------------------------------------------------------------------------
// Structural transformations
// ---------------------------------------------------------------------------

matrix matrix::transpose() const {
    matrix m(col, row);
    for (int i = 0; i < col; i++)
        for (int j = 0; j < row; j++)
            m.arr[i * row + j] = arr[j * col + i];
    return m;
}

matrix matrix::adjoint() {
    if (row != col)
        throw std::runtime_error("Adjoint is defined only for square matrix.");
    matrix m(*this);
    matrix adj(row, col);
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            adj.arr[i * col + j] = m.cofactor(i, j);
    fpg(adj);
    return adj.transpose();
}

matrix matrix::inverse() const {
    if (row != col)
        throw std::runtime_error("Inverse is defined only for square matrix.");
    const int m_sz = row;
    matrix mat(*this);
    matrix inv = identity(m_sz);

    for (int k = 0; k < m_sz; k++) {
        // Partial pivot
        int z = k;
        while (z < m_sz && std::abs(mat.arr[z * m_sz + k]) <= EPS) z++;
        if (z == m_sz)
            throw std::runtime_error("Inverse is defined only for Non-Singular matrix.");
        if (z != k) {
            mat = mat.row_swap(k, z);
            inv = inv.row_swap(k, z);
        }
        const double divisor = mat.arr[k * m_sz + k];
        if (std::abs(divisor) <= EPS)
            throw std::runtime_error("Inverse is defined only for Non-Singular matrix.");
        for (int j = 0; j < m_sz; j++) {
            mat.arr[k * m_sz + j] /= divisor;
            inv.arr[k * m_sz + j] /= divisor;
        }
        for (int i = 0; i < m_sz; i++) {
            if (i == k) continue;
            const double multiplier = mat.arr[i * m_sz + k];
            if (std::abs(multiplier) < 1e-15) continue;
            for (int j = 0; j < m_sz; j++) {
                mat.arr[i * m_sz + j] -= multiplier * mat.arr[k * m_sz + j];
                inv.arr[i * m_sz + j] -= multiplier * inv.arr[k * m_sz + j];
            }
        }
    }
    neg_zero(inv);
    fpg(inv);
    return inv;
}

matrix matrix::solve() {
    if (col - row != 1)
        throw std::runtime_error("The dimension of Augmented Matrix must be (N x N+1).");
    matrix mat(*this);
    mat.gaussian();
    matrix solution(row, 1);
    for (int i = row - 1; i >= 0; i--) {
        double subtractor = 0;
        for (int j = i + 1; j < row; j++)
            subtractor += solution.arr[j] * mat.arr[i * col + j];
        solution.arr[i] = mat.arr[i * col + row] - subtractor;
    }
    fpg(solution);
    return solution;
}

// ---------------------------------------------------------------------------
// Gram-Schmidt / QR
// ---------------------------------------------------------------------------

matrix matrix::orthogonalize() {
    matrix Q = *this;
    for (int j = 0; j < col; j++) {
        for (int i = 0; i < j; i++) {
            double dot = 0.0, norm_sq = 0.0;
            for (int k = 0; k < row; k++) {
                dot     += Q.arr[k * col + j] * Q.arr[k * col + i];
                norm_sq += Q.arr[k * col + i] * Q.arr[k * col + i];
            }
            const double proj = (norm_sq > EPS) ? dot / norm_sq : 0.0;
            for (int k = 0; k < row; k++)
                Q.arr[k * col + j] -= proj * Q.arr[k * col + i];
        }
        double norm_sq_j = 0.0;
        for (int k = 0; k < row; k++)
            norm_sq_j += Q.arr[k * col + j] * Q.arr[k * col + j];
        if (norm_sq_j <= EPS)
            for (int k = 0; k < row; k++)
                Q.arr[k * col + j] = 0.0;
    }
    fpg(Q);
    return Q;
}

matrix matrix::orthonormalize() {
    matrix Q = this->orthogonalize();
    for (int j = 0; j < col; j++) {
        double norm = 0.0;
        for (int k = 0; k < row; k++)
            norm += Q.arr[k * col + j] * Q.arr[k * col + j];
        norm = std::sqrt(norm);
        if (norm > EPS)
            for (int k = 0; k < row; k++)
                Q.arr[k * col + j] /= norm;
        else
            for (int k = 0; k < row; k++)
                Q.arr[k * col + j] = 0.0;
    }
    fpg(Q);
    return Q;
}

matrix matrix::qr_decomp_q() { return orthonormalize(); }

matrix matrix::qr_decomp_r() {
    matrix Q = qr_decomp_q();
    matrix R(col, col);
    for (int j = 0; j < col; j++)
        for (int i = 0; i <= j; i++) {
            double val = 0.0;
            for (int k = 0; k < row; k++)
                val += Q.arr[k * col + i] * arr[k * col + j];
            R.arr[i * col + j] = val;
        }
    fpg(R);
    return R;
}

// ---------------------------------------------------------------------------
// LU decomposition
// ---------------------------------------------------------------------------

std::tuple<matrix, matrix, matrix> matrix::lu_decomposition() const {
    if (row != col)
        throw std::runtime_error("LU decomposition requires a square matrix.");
    const size_t n = static_cast<size_t>(row);
    matrix P(row, col);
    matrix L(row, col);
    matrix U = *this;

    for (size_t i = 0; i < n; i++) {
        P.arr[i * n + i] = 1.0;
        L.arr[i * n + i] = 1.0;
    }
    for (size_t i = 0; i < n; i++) {
        size_t pivot_row = i;
        double max_val = std::abs(U.arr[i * n + i]);
        for (size_t k = i + 1; k < n; k++)
            if (std::abs(U.arr[k * n + i]) > max_val) {
                max_val = std::abs(U.arr[k * n + i]);
                pivot_row = k;
            }
        if (max_val < 1e-9)
            throw std::runtime_error("Matrix is singular; LU decomposition failed.");
        if (pivot_row != i) {
            for (size_t j = 0; j < n; j++) {
                std::swap(U.arr[i * n + j], U.arr[pivot_row * n + j]);
                std::swap(P.arr[i * n + j], P.arr[pivot_row * n + j]);
            }
            for (size_t j = 0; j < i; j++)
                std::swap(L.arr[i * n + j], L.arr[pivot_row * n + j]);
        }
        for (size_t k = i + 1; k < n; k++) {
            double factor = U.arr[k * n + i] / U.arr[i * n + i];
            L.arr[k * n + i] = factor;
            for (size_t j = i; j < n; j++)
                U.arr[k * n + j] -= factor * U.arr[i * n + j];
        }
    }
    neg_zero(L); fpg(L);
    neg_zero(U); fpg(U);
    return {P, L, U};
}

// ---------------------------------------------------------------------------
// Triangular checks (tolerant)
// ---------------------------------------------------------------------------

bool matrix::check_upper_tri() {
    for (int j = 0; j < col; j++)
        for (int i = j + 1; i < row; i++)
            // Tolerant zero-check: a value is "not zero" only if |v| > APPROX_TOL
            if (!approx_equal(arr[i * col + j], 0.0, APPROX_TOL)) return false;
    return true;
}

bool matrix::check_lower_tri() {
    for (int i = 0; i < row; i++)
        for (int j = i + 1; j < col; j++)
            if (!approx_equal(arr[i * col + j], 0.0, APPROX_TOL)) return false;
    return true;
}

// ---------------------------------------------------------------------------
// Eigenvalues (unshifted QR iteration)
// ---------------------------------------------------------------------------

matrix matrix::char_poly() {
    if (row != col)
        throw std::runtime_error("Characteristic Polynomial is defined only for square matrix.");
    int n = row;
    matrix poly(1, n + 1);
    matrix id = identity(n);
    matrix b_k_1 = id;
    poly.arr[0] = 1.0;
    for (int i = 1; i <= n; i++) {
        matrix a = (*this) * b_k_1;
        double c_k = -a.trace() / i;
        matrix b_k = a + id * c_k;
        poly.arr[i] = c_k;
        b_k_1 = b_k;
    }
    return poly;
}

matrix matrix::eigenvalues() {
    if (row != col)
        throw std::runtime_error("Eigenvalues are defined only for square matrix.");
    matrix m(*this);
    const int max_iter = 1000;
    bool converged = false;
    for (int iter = 0; iter < max_iter; iter++) {
        matrix q = m.qr_decomp_q();
        matrix r = m.qr_decomp_r();
        m = r * q;
        fpg(m);
        converged = true;
        for (int j = 0; j < m.col; j++)
            for (int i = j + 1; i < m.row; i++)
                if (std::abs(m.arr[i * m.col + j]) > EPS) { converged = false; break; }
        if (!converged) continue;
        break;
    }
    if (!converged)
        throw std::runtime_error(
            "QR iteration failed to converge after 1000 iterations. "
            "The matrix may have complex eigenvalues or require shift techniques.");
    matrix eigen(m.row, 1);
    for (int i = 0; i < m.row; i++)
        eigen.arr[i] = m.arr[i * m.col + i];
    return eigen;
}

double matrix::norm() const {
    return std::sqrt(inner_product(*this, *this));
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

matrix matpow(matrix mat, long long expo) {
    if (mat.row != mat.col)
        throw std::runtime_error("Matrix power is defined only for square matrices.");
    if (expo < 0)
        throw std::runtime_error("Matrix power is defined only for non-negative exponents.");
    matrix res = identity(mat.row);
    while (expo > 0) {
        if (expo & 1) res = res * mat;
        mat = mat * mat;
        expo >>= 1;
    }
    return res;
}

bool check_ortho(const matrix& mat, double tol) {
    if (mat.row != mat.col) return false; // non-square cannot be orthogonal
    try {
        // Q is orthogonal iff Q^T == Q^-1, checked element-wise within tol.
        return matrices_approx_equal(mat.transpose(), mat.inverse(), tol);
    } catch (const std::runtime_error&) {
        // Singular matrix has no inverse, so cannot be orthogonal.
        return false;
    }
}

bool check_unitary(const matrix& mat, double tol) {
    if (mat.row != mat.col) return false;
    try {
        matrix tmp(mat);
        return matrices_approx_equal(mat.transpose(), tmp.adjoint(), tol);
    } catch (const std::runtime_error&) {
        return false;
    }
}

double inner_product(const matrix& a, const matrix& b) {
    if (a.row != b.row || a.col != b.col)
        throw std::runtime_error("Vectors don't belong same vector space.");
    double sum = 0.0;
    const size_t n = static_cast<size_t>(a.row) * static_cast<size_t>(a.col);
    for (size_t i = 0; i < n; i++)
        sum += a.arr[i] * b.arr[i];
    return sum;
}

double angle(const matrix& a, const matrix& b) {
    double norm_a = a.norm();
    double norm_b = b.norm();
    if (norm_a < 1e-9 || norm_b < 1e-9)
        throw std::runtime_error(
            "Cannot compute vector angle involving a zero-norm matrix space.");
    double cos_theta = inner_product(a, b) / (norm_a * norm_b);
    // Guard against minor floating-point overshoots (e.g. 1.0000000002)
    if (cos_theta >  1.0) cos_theta =  1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    return std::acos(cos_theta);
}

} // namespace linalg
