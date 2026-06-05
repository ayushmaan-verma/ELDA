#include "elda/matrix.hpp"
#include <algorithm>
#include <tuple>
#include <stdexcept>
#include <cmath>
#include <ostream>

namespace linalg {

void matrix::input() {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            std::cin >> arr[i * col + j];
        }
    }
}

void matrix::print() {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            std::cout << arr[i * col + j] << " ";
        }
        std::cout << std::endl;
    }
}

std::ostream& operator<<(std::ostream& os, const matrix& m) {
    for (int i = 0; i < m.get_rows(); i++) {
        for (int j = 0; j < m.get_cols(); j++) {
            os << m(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}

matrix identity(int n) {
    if (n < 0) {
        throw std::runtime_error("Identity matrix size must be non-negative.");
    }
    matrix m(n, n);
    for (int i = 0; i < n; i++) {
        m.arr[i * n + i] = 1.0;
    }
    return m;
}

void neg_zero(matrix& m) {
    size_t total_elements = m.get_rows() * m.get_cols();
    for (size_t i = 0; i < total_elements; i++) {
        if (m.arr[i] == -0.0) {
            m.arr[i] = 0.0;
        }
    }
}

void fpg(matrix& m) {
    size_t total_elements = m.get_rows() * m.get_cols();
    for (size_t i = 0; i < total_elements; i++) {
        if (std::abs(m.arr[i]) <= EPS) {
            m.arr[i] = 0.0;
        }
    }
}

matrix matrix::operator=(const matrix& m2) {
    if ((row != m2.row) || (col != m2.col)) {
        throw std::runtime_error("Assignment Operator : Dimension Mismatch in LHS & RHS !!!");
    }
    arr = m2.arr;
    fpg(*this);
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
    fpg(m3);
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
    fpg(m3);
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

bool operator==(const matrix& m1, const matrix& m2) {
    return m1.arr == m2.arr;
}

bool shape_comp(const matrix& m1, const matrix& m2) {
    return (m1.row == m2.row) && (m1.col == m2.col);
}

matrix matrix::get_row(int r) {
    matrix rw(row, col);
    for (int j = 0; j < col; j++) {
        rw.arr[r * col + j] = arr[r * col + j];
    }
    return rw;
}

matrix matrix::get_col(int c) {
    matrix cn(row, col);
    for (int i = 0; i < row; i++) {
        cn.arr[i * col + c] = arr[i * col + c];
    }
    return cn;
}

matrix matrix::get_row_vec(int r) {
    matrix rw(1, col);
    for (int j = 0; j < col; j++) {
        rw.arr[j] = arr[r * col + j];
    }
    return rw;
}

matrix matrix::get_col_vec(int c) {
    matrix cn(row, 1);
    for (int i = 0; i < row; i++) {
        cn.arr[i] = arr[i * col + c];
    }
    return cn;
}

void matrix::replace_row(int r, const matrix& rw) {
    for (int j = 0; j < col; j++) {
        arr[r * col + j] = rw.arr[j];
    }
}

void matrix::replace_col(int c, const matrix& cn) {
    for (int i = 0; i < row; i++) {
        arr[i * col + c] = cn.arr[i];
    }
}

matrix matrix::row_op(int i, double coeff, int j) {
    matrix m(*this);
    for (int k = 0; k < col; k++) {
        m.arr[i * col + k] += coeff * arr[j * col + k];
    }
    fpg(m);
    return m;
}

matrix matrix::col_op(int i, double coeff, int j) {
    matrix m(*this);
    for (int k = 0; k < row; k++) {
        m.arr[k * col + i] += coeff * m.arr[k * col + j];
    }
    fpg(m);
    return m;
}

matrix matrix::row_swap(int i, int j) {
    matrix m(*this);
    for (int k = 0; k < col; k++) {
        std::swap(m.arr[i * col + k], m.arr[j * col + k]);
    }
    fpg(m);
    return m;
}

matrix matrix::col_swap(int i, int j) {
    matrix m(*this);
    for (int k = 0; k < row; k++) {
        std::swap(m.arr[k * col + i], m.arr[k * col + j]);
    }
    fpg(m);
    return m;
}

matrix matrix::row_multi(int i, double factor) {
    matrix m(*this);
    for (int k = 0; k < col; k++) {
        m.arr[i * col + k] *= factor;
    }
    fpg(m);
    return m;
}

matrix matrix::col_multi(int i, double factor) {
    matrix m(*this);
    for (int k = 0; k < row; k++) {
        m.arr[k * col + i] *= factor;
    }
    fpg(m);
    return m;
}

int matrix::echelon() {
    int swaps = 1;
    for (int k = 0; k < std::min(row, col); k++) {
        int z = k;
        while (std::abs(arr[z * col + k]) < 1e-9) {
            if (z + 1 == std::min(row, col)) {
                break;
            }
            z++;
        }
        if (k != z) {
            *this = row_swap(k, z);
            swaps *= -1;
        }
        if (std::abs(arr[k * col + k]) < 1e-9) {
            continue;
        }
        for (int i = k + 1; i < row; i++) {
            const double multiplier = arr[i * col + k] / arr[k * col + k];
            *this = row_op(i, -multiplier, k);
        }
    }
    neg_zero(*this);
    fpg(*this);
    return swaps;
}

int matrix::gaussian() {
    int swaps = 1;
    for (int k = 0; k < std::min(row, col); k++) {
        int z = k;
        while (std::abs(arr[z * col + k]) < 1e-9) {
            if (z + 1 == std::min(row, col)) {
                break;
            }
            z++;
        }
        if (k != z) {
            *this = row_swap(k, z);
            swaps *= -1;
        }
        const double divisor = arr[k * col + k];
        if (std::abs(divisor) < 1e-9) {
            continue;
        }
        *this = this->row_multi(k, 1.0 / divisor);
        for (int i = k + 1; i < row; i++) {
            const double multiplier = arr[i * col + k];
            *this = row_op(i, -multiplier, k);
        }
    }
    neg_zero(*this);
    fpg(*this);
    return swaps;
}

int matrix::gauss_jordan() {
    const int swaps = this->gaussian();
    for (int k = std::min(row, col) - 1; k > 0; k--) {
        for (int i = k - 1; i >= 0; i--) {
            const double multiplier = arr[i * col + k];
            *this = row_op(i, -multiplier, k);
        }
    }
    fpg(*this);
    return swaps;
}

int matrix::canonical() {
    const int swaps = this->gaussian();
    for (int k = 0; k < std::min(row, col); k++) {
        for (int j = k + 1; j < col; j++) {
            const double multiplier = arr[k * col + j];
            *this = col_op(j, -multiplier, k);
        }
    }
    fpg(*this);
    return swaps;
}

int matrix::rank() {
    int rank_cnt = 0;
    matrix m(*this);
    m.gaussian();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (m.arr[i * col + j] != 0.0) {
                rank_cnt++;
                break;
            }
        }
    }
    return rank_cnt;
}

int matrix::free_variable() {
    if (col != row) {
        throw std::runtime_error("The dimension of Coefficient Matrix must be N x N.");
    }
    return row - this->rank();
}

double matrix::trace() {
    if (row != col) {
        throw std::runtime_error("Trace is defined only for square matrix.");
    }
    double sum = 0;
    for (int i = 0; i < row; i++) {
        sum += arr[i * col + i];
    }
    return sum;
}

double matrix::diag_product() {
    double prod = 1.0;
    for (int i = 0; i < std::min(row, col); i++) {
        prod *= this->arr[i * col + i];
    }
    return prod;
}

double matrix::det() {
    if (row != col) {
        throw std::runtime_error("Determinant is defined only for square matrix.");
    }
    matrix m(*this);
    double determinant = m.echelon();
    determinant *= m.diag_product();
    if (determinant == -0.0) {
        determinant = 0.0;
    }
    return determinant;
}

double matrix::minor(int x, int y) {
    if (row != col) {
        throw std::runtime_error("Cofactor & Minor defined only for square matrix.");
    }
    matrix mat(row - 1, col - 1);
    for (int i = 0, r = 0; i < row; i++) {
        if (i == x) continue;
        for (int j = 0, c = 0; j < col; j++) {
            if (j == y) continue;
            mat.arr[r * (col - 1) + c] = arr[i * col + j];
            c++;
        }
        r++;
    }
    return mat.det();
}

double matrix::cofactor(int x, int y) {
    int check = x + y + 2;
    if (check % 2) {
        return -this->minor(x, y);
    }
    return this->minor(x, y);
}

matrix matrix::transpose() {
    matrix m(col, row);
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < row; j++) {
            m.arr[i * row + j] = arr[j * col + i];
        }
    }
    return m;
}

matrix matrix::adjoint() {
    if (row != col) {
        throw std::runtime_error("Adjoint is defined only for square matrix.");
    }
    matrix m(*this);
    matrix adj(row, col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            adj.arr[i * col + j] = m.cofactor(i, j);
        }
    }
    fpg(adj);
    return adj.transpose();
}

matrix matrix::inverse() {
    if (row != col) {
        throw std::runtime_error("Inverse is defined only for square matrix.");
    }
    const int m = row;
    matrix mat(*this);
    matrix inv = identity(m);

    for (int k = 0; k < m; k++) {
        int z = k;
        while (std::abs(mat.arr[z * m + k]) < 1e-9) {
            if (z + 1 == m) break;
            z++;
        }
        if (z != k) {
            mat = mat.row_swap(k, z);
            inv = inv.row_swap(k, z);
        }
        const double divisor = mat.arr[k * m + k];
        mat = mat.row_multi(k, 1.0 / divisor);
        inv = inv.row_multi(k, 1.0 / divisor);

        for (int i = k + 1; i < m; i++) {
            const double multiplier = mat.arr[i * m + k];
            mat = mat.row_op(i, -multiplier, k);
            inv = inv.row_op(i, -multiplier, k);
        }
    }
    for (int k = m - 1; k > 0; k--) {
        for (int i = k - 1; i >= 0; i--) {
            const double multiplier = mat.arr[i * m + k];
            mat = mat.row_op(i, -multiplier, k);
            inv = inv.row_op(i, -multiplier, k);
        }
    }
    for (int i = 0; i < m; i++) {
        if (std::abs(mat.arr[i * m + i] - 1.0) > 1e-5) {
            throw std::runtime_error("Inverse is defined only for Non-Singular matrix.");
        }
    }
    neg_zero(inv);
    fpg(inv);
    return inv;
}

matrix matrix::solve() {
    if (col - row != 1) {
        throw std::runtime_error("The dimension of Augmented Matrix must be (N x N+1).");
    }
    matrix mat(*this);
    mat.gaussian();
    matrix solution(row, 1);
    for (int i = row - 1; i >= 0; i--) {
        double subtractor = 0;
        for (int j = i + 1; j < row; j++) {
            subtractor += solution.arr[j] * mat.arr[i * col + j];
        }
        solution.arr[i] = mat.arr[i * col + row] - subtractor;
    }
    fpg(solution);
    return solution;
}

matrix matrix::orthogonalize() {
    matrix Q = *this;
    for (int j = 0; j < col; j++) {
        for (int i = 0; i < j; i++) {
            double dot_product = 0.0;
            double norm_sq_i = 0.0;
            for (int k = 0; k < row; k++) {
                dot_product += Q.arr[k * col + j] * Q.arr[k * col + i];
                norm_sq_i += Q.arr[k * col + i] * Q.arr[k * col + i];
            }
            double projection_coeff = 0.0;
            if (norm_sq_i > 1e-9) {
                projection_coeff = dot_product / norm_sq_i;
            }
            for (int k = 0; k < row; k++) {
                Q.arr[k * col + j] -= projection_coeff * Q.arr[k * col + i];
            }
        }
    }
    fpg(Q);
    return Q;
}

matrix matrix::orthonormalize() {
    matrix Q = this->orthogonalize();
    for (int j = 0; j < col; j++) {
        double norm_val = 0.0;
        for (int k = 0; k < row; k++) {
            norm_val += Q.arr[k * col + j] * Q.arr[k * col + j];
        }
        norm_val = std::sqrt(norm_val);
        if (norm_val > 1e-9) {
            for (int k = 0; k < row; k++) {
                Q.arr[k * col + j] /= norm_val;
            }
        } else {
            for (int k = 0; k < row; k++) {
                Q.arr[k * col + j] = 0.0;
            }
        }
    }
    fpg(Q);
    return Q;
}

matrix matrix::qr_decomp_q() {
    return orthonormalize();
}

matrix matrix::qr_decomp_r() {
    matrix Q = qr_decomp_q();
    matrix R(col, col);
    for (int j = 0; j < col; j++) {
        for (int i = 0; i <= j; i++) {
            double val = 0.0;
            for (int k = 0; k < row; k++) {
                val += Q.arr[k * col + i] * arr[k * col + j];
            }
            R.arr[i * col + j] = val;
        }
    }
    fpg(R);
    return R;
}

std::tuple<matrix, matrix, matrix> matrix::lu_decomposition() const {
    if (row != col) {
        throw std::runtime_error("LU decomposition requires a square matrix.");
    }
    size_t n = static_cast<size_t>(row);
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
        for (size_t k = i + 1; k < n; k++) {
            if (std::abs(U.arr[k * n + i]) > max_val) {
                max_val = std::abs(U.arr[k * n + i]);
                pivot_row = k;
            }
        }
        if (max_val < 1e-9) {
            throw std::runtime_error("Matrix is singular; LU decomposition failed.");
        }
        if (pivot_row != i) {
            for (size_t j = 0; j < n; j++) {
                std::swap(U.arr[i * n + j], U.arr[pivot_row * n + j]);
                std::swap(P.arr[i * n + j], P.arr[pivot_row * n + j]);
            }
            for (size_t j = 0; j < i; j++) {
                std::swap(L.arr[i * n + j], L.arr[pivot_row * n + j]);
            }
        }
        for (size_t k = i + 1; k < n; k++) {
            double factor = U.arr[k * n + i] / U.arr[i * n + i];
            L.arr[k * n + i] = factor;
            for (size_t j = i; j < n; j++) {
                U.arr[k * n + j] -= factor * U.arr[i * n + j];
            }
        }
    }
    neg_zero(L); fpg(L);
    neg_zero(U); fpg(U);
    return {P, L, U};
}

bool matrix::check_upper_tri() {
    for (int j = 0; j < col; j++) {
        for (int i = j + 1; i < row; i++) {
            if (arr[i * col + j] != 0.0) return false;
        }
    }
    return true;
}

bool matrix::check_lower_tri() {
    for (int i = 0; i < row; i++) {
        for (int j = i + 1; j < col; j++) {
            if (arr[i * col + j] != 0.0) return false;
        }
    }
    return true;
}

matrix matrix::char_poly() {
    if (this->row != this->col) {
        throw std::runtime_error("Characteristic Polynomial is defined only for square matrix.");
    }
    int n = this->row;
    matrix poly(1, n + 1);
    matrix b_k(n, n);
    matrix id = identity(n);
    matrix b_k_1 = id;
    poly.arr[0] = 1.0;
    for (int i = 1; i <= n; i++) {
        matrix a = (*this) * b_k_1;
        double c_k = -a.trace() / i;
        b_k = a + id * c_k;
        poly.arr[i] = c_k;
        b_k_1 = b_k;
    }
    return poly;
}

matrix matrix::eigenvalues() {
    if (row != col) {
        throw std::runtime_error("Eigenvalues are defined only for square matrix.");
    }
    matrix m(*this);
    int max_iter = 1000;
    int iter = 0;
    bool converged = false;

    while (iter < max_iter) {
        matrix q = m.qr_decomp_q();
        matrix r = m.qr_decomp_r();
        m = r * q;
        fpg(m);

        converged = true;
        for (int j = 0; j < m.col; j++) {
            for (int i = j + 1; i < m.row; i++) {
                if (std::abs(m.arr[i * m.col + j]) > EPS) {
                    converged = false;
                    break;
                }
            }
            if (!converged) break;
        }
        if (converged) break;
        iter++;
    }
    if (!converged) {
        throw std::runtime_error("QR iteration failed to converge.");
    }
    matrix eigen(row, 1);
    for (int i = 0; i < row; i++) {
        eigen.arr[i] = m.arr[i * col + i];
    }
    return eigen;
}

double matrix::norm() {
    return sqrt(inner_product(*this, *this));
}

matrix matpow(matrix mat, long long expo) {
    if (mat.row != mat.col) {
        throw std::runtime_error("Matrix exponentiation is defined only for square matrix.");
    }
    if (expo < 0) {
        throw std::runtime_error("Matrix exponentiation requires a non-negative exponent.");
    }
    matrix res = identity(mat.row);
    while (expo > 0) {
        if (expo & 1) res = res * mat;
        mat = mat * mat;
        expo >>= 1;
    }
    return res;
}

bool check_ortho(const matrix& mat) {
    matrix m(mat);
    matrix mt = m.transpose();
    matrix inv = m.inverse();
    const size_t total_elements = mt.get_rows() * mt.get_cols();
    for (size_t i = 0; i < total_elements; i++) {
        if (std::abs(mt.arr[i] - inv.arr[i]) > EPS) {
            return false;
        }
    }
    return true;
}

bool check_unitary(const matrix& mat) {
    matrix m(mat);
    matrix mt = m.transpose();
    matrix adj = m.adjoint();
    const size_t total_elements = mt.get_rows() * mt.get_cols();
    for (size_t i = 0; i < total_elements; i++) {
        if (std::abs(mt.arr[i] - adj.arr[i]) > EPS) {
            return false;
        }
    }
    return true;
}

double inner_product(const matrix& a, const matrix& b) {
    if ((a.row != b.row) || (a.col != b.col)) {
        throw std::runtime_error("Vectors don't belong same vector space.");
    }
    double sum = 0.0;
    const size_t total_elements = a.get_rows() * a.get_cols();
    for (size_t i = 0; i < total_elements; i++) {
        sum += a.arr[i] * b.arr[i];
    }
    return sum;
}

double angle(const matrix& a, const matrix& b) {
    double l_a = std::sqrt(inner_product(a, a));
    double l_b = std::sqrt(inner_product(b, b));
    if (l_a <= EPS || l_b <= EPS) {
        throw std::runtime_error("Angle undefined for zero-length vectors.");
    }
    double cosine = inner_product(a, b) / (l_a * l_b);
    cosine = std::max(-1.0, std::min(1.0, cosine));
    return acos(cosine);
}

} // namespace linalg