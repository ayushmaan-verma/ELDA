#include "elda/matrix.hpp"

#include <algorithm>

namespace linalg {

void matrix::input() {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            std::cin >> arr[i][j];
        }
    }
}

void matrix::print() {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

matrix identity(int n) {
    matrix m(n, n);
    for (int i = 0; i < n; i++) {
        m.arr[i][i] = 1;
    }
    return m;
}

void neg_zero(matrix& m) {
    for (int i = 0; i < m.row; i++) {
        for (int j = 0; j < m.col; j++) {
            if (m.arr[i][j] == -0) {
                m.arr[i][j] = 0;
            }
        }
    }
}

void fpg(matrix& m) {
    // Snap tiny values to zero to avoid floating-point residue.
    for (int i = 0; i < m.row; i++) {
        for (int j = 0; j < m.col; j++) {
            if (std::abs(m.get_element(i, j)) <= EPS) {
                *m.ref_element(i, j) = 0;
            }
        }
    }
}

matrix matrix::operator=(matrix m2) {
    if ((row != m2.row) || (col != m2.col)) {
        throw std::runtime_error("Assignment Operator : Dimension Mismatch in LHS & RHS !!!");
    }
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            arr[i][j] = m2.arr[i][j];
        }
    }
    fpg(*this);
    return *this;
}

matrix matrix::operator+(matrix m2) {
    if ((row != m2.row) || (col != m2.col)) {
        throw std::runtime_error("Addition of two matrices is defined only if both of them are of same dimensions.");
    }
    matrix m3(row, col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            m3.arr[i][j] = this->arr[i][j] + m2.arr[i][j];
        }
    }
    fpg(m3);
    return m3;
}

matrix matrix::operator-(matrix m2) {
    if ((row != m2.row) || (col != m2.col)) {
        throw std::runtime_error("Subtraction of two matrices is defined only if both of them are of same dimensions.");
    }
    matrix m3(row, col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            m3.arr[i][j] = this->arr[i][j] - m2.arr[i][j];
        }
    }
    fpg(m3);
    return m3;
}

matrix matrix::operator*(matrix m2) {
    if (col != m2.row) {
        throw std::runtime_error(
            "Matrix multiplication can be performed only if number of columns of 1st matrix matches the number of rows of 2nd matrix.");
    }
    matrix m3(row, m2.col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < m2.col; j++) {
            for (int k = 0; k < col; k++) {
                m3.arr[i][j] += this->arr[i][k] * m2.arr[k][j];
            }
        }
    }
    fpg(m3);
    return m3;
}

matrix matrix::operator*(double d) {
    matrix m3(*this);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            m3.arr[i][j] *= d;
        }
    }
    fpg(m3);
    return m3;
}

bool operator==(matrix m1, matrix m2) {
    return m1.arr == m2.arr;
}

bool shape_comp(matrix m1, matrix m2) {
    return (m1.row == m2.row) && (m1.col == m2.col);
}

matrix matrix::get_row(int r) {
    matrix rw(row, col);
    rw.arr[r] = arr[r];
    return rw;
}

matrix matrix::get_col(int c) {
    matrix cn(row, col);
    for (int i = 0; i < row; i++) {
        cn.arr[i][c] = arr[i][c];
    }
    return cn;
}

matrix matrix::get_row_vec(int r) {
    matrix rw(1, col);
    rw.arr[0] = arr[r];
    return rw;
}

matrix matrix::get_col_vec(int c) {
    matrix cn(row, 1);
    for (int i = 0; i < row; i++) {
        cn.arr[i][0] = arr[i][c];
    }
    return cn;
}

void matrix::replace_row(int r, matrix rw) {
    arr[r] = rw.arr[0];
}

void matrix::replace_col(int c, matrix cn) {
    for (int i = 0; i < row; i++) {
        arr[i][c] = cn.arr[i][0];
    }
}

matrix matrix::row_op(int i, double coeff, int j) {
    matrix m(*this);
    for (int k = 0; k < col; k++) {
        m.arr[i][k] += coeff * arr[j][k];
    }
    fpg(m);
    return m;
}

matrix matrix::col_op(int i, double coeff, int j) {
    matrix m(*this);
    for (int k = 0; k < row; k++) {
        m.arr[k][i] += coeff * m.arr[k][j];
    }
    fpg(m);
    return m;
}

matrix matrix::row_swap(int i, int j) {
    matrix m(*this);
    std::swap(m.arr[i], m.arr[j]);
    fpg(m);
    return m;
}

matrix matrix::col_swap(int i, int j) {
    matrix m(*this);
    for (int k = 0; k < row; k++) {
        std::swap(m.arr[k][i], m.arr[k][j]);
    }
    fpg(m);
    return m;
}

matrix matrix::row_multi(int i, double factor) {
    matrix m(*this);
    for (int k = 0; k < col; k++) {
        m.arr[i][k] *= factor;
    }
    fpg(m);
    return m;
}

matrix matrix::col_multi(int i, double factor) {
    matrix m(*this);
    for (int k = 0; k < row; k++) {
        m.arr[k][i] *= factor;
    }
    fpg(m);
    return m;
}

int matrix::echelon() {
    int swaps = 1;
    for (int k = 0; k < std::min(row, col); k++) {
        int z = k;
        while (!static_cast<bool>(arr[z][k])) {
            if (z + 1 == std::min(row, col)) {
                break;
            }
            z++;
        }
        if (k != z) {
            *this = row_swap(k, z);
            swaps *= -1;
        }
        if (arr[k][k] == 0) {
            continue;
        }
        // Eliminate entries below the current pivot.
        for (int i = k + 1; i < row; i++) {
            const double multiplier = arr[i][k] / arr[k][k];
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
        while (!static_cast<bool>(arr[z][k])) {
            if (z + 1 == std::min(row, col)) {
                break;
            }
            z++;
        }
        if (k != z) {
            *this = row_swap(k, z);
            swaps *= -1;
        }
        const double divisor = arr[k][k];
        if (divisor == 0) {
            continue;
        }
        *this = this->row_multi(k, 1 / divisor);
        // Normalize the pivot row, then eliminate entries below it.
        for (int i = k + 1; i < row; i++) {
            const double multiplier = arr[i][k];
            *this = row_op(i, -multiplier, k);
        }
    }
    neg_zero(*this);
    fpg(*this);
    return swaps;
}

int matrix::gauss_jordan() {
    const int swaps = this->gaussian();
    // Clear the entries above each pivot after forward elimination.
    for (int k = std::min(row, col) - 1; k > 0; k--) {
        for (int i = k - 1; i >= 0; i--) {
            const double multiplier = arr[i][k];
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
            const double multiplier = arr[k][j];
            *this = col_op(j, -multiplier, k);
        }
    }
    fpg(*this);
    return swaps;
}

int matrix::rank() {
    int rank = 0;
    matrix m(*this);
    m.gaussian();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (m.arr[i][j] != 0) {
                rank++;
                break;
            }
        }
    }
    return rank;
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
    for (int i = 0; i < std::min(row, col); i++) {
        sum += arr[i][i];
    }
    return sum;
}

double matrix::diag_product() {
    double prod = 1;
    for (int i = 0; i < std::min(row, col); i++) {
        prod *= this->arr[i][i];
    }
    return prod;
}

double matrix::det() {
    if (row != col) {
        throw std::runtime_error("Determinant is defined only for square matrix.");
    }
    matrix m(*this);
    double det = m.echelon();
    det *= m.diag_product();
    if (det == -0) {
        det = 0;
    }
    return det;
}

double matrix::minor(int x, int y) {
    if (row != col) {
        throw std::runtime_error("Cofactor & Minor of an element is defined only for square matrix.");
    }
    matrix mat(row - 1, col - 1);
    for (int i = 0, r = 0; i < row; i++) {
        if (i == x) {
            continue;
        }
        for (int j = 0, c = 0; j < col; j++) {
            if (j == y) {
                continue;
            }
            mat.arr[r][c] = arr[i][j];
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
            m.arr[i][j] = arr[j][i];
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
            adj.arr[i][j] = m.cofactor(i, j);
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

    // Apply Gauss-Jordan elimination to [mat | inv].
    for (int k = 0; k < m; k++) {
        int z = k;
        while (!static_cast<bool>(mat.arr[z][k])) {
            if (z + 1 == m) {
                break;
            }
            z++;
        }

        if (z != k) {
            mat = mat.row_swap(k, z);
            inv = inv.row_swap(k, z);
        }

        const double divisor = mat.arr[k][k];
        mat = mat.row_multi(k, 1 / divisor);
        inv = inv.row_multi(k, 1 / divisor);

        for (int i = k + 1; i < m; i++) {
            const double multiplier = mat.arr[i][k];
            mat = mat.row_op(i, -multiplier, k);
            inv = inv.row_op(i, -multiplier, k);
        }
    }

    for (int k = m - 1; k > 0; k--) {
        for (int i = k - 1; i >= 0; i--) {
            const double multiplier = mat.arr[i][k];
            mat = mat.row_op(i, -multiplier, k);
            inv = inv.row_op(i, -multiplier, k);
        }
    }

    neg_zero(inv);

    for (int i = 0; i < m; i++) {
        if (mat.arr[i][i] != 1) {
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
    // Perform back-substitution on the upper-triangular system.
    for (int i = row - 1; i >= 0; i--) {
        double subtractor = 0;
        for (int j = i + 1; j < row; j++) {
            subtractor += solution.arr[j][0] * mat.arr[i][j];
        }
        solution.arr[i][0] = mat.arr[i][row] - subtractor;
        subtractor = 0;
    }
    fpg(solution);
    return solution;
}

matrix matrix::orthogonalize() {
    matrix ot(row, col);
    matrix v = get_col_vec(0);
    matrix w = v;
    for (int i = 0; i < col; i++) {
        ot.replace_col(i, w);
        v = get_col_vec(i + 1);
        for (int j = 0; j < i + 1; j++) {
            v = v - ot.get_col_vec(j) *
                        (inner_product(v, ot.get_col_vec(j)) / inner_product(ot.get_col_vec(j), ot.get_col_vec(j)));
        }
        w = v;
    }
    return ot;
}

matrix matrix::orthonormalize() {
    matrix ot(row, col);
    matrix v = get_col_vec(0);
    matrix w = v * (1 / v.norm());
    for (int i = 0; i < col; i++) {
        ot.replace_col(i, w);
        v = get_col_vec(i + 1);
        for (int j = 0; j < i + 1; j++) {
            v = v - ot.get_col_vec(j) *
                        (inner_product(v, ot.get_col_vec(j)) / inner_product(ot.get_col_vec(j), ot.get_col_vec(j)));
        }
        w = v * (1 / v.norm());
    }
    return ot;
}

matrix matrix::qr_decomp_q() {
    return orthonormalize();
}

matrix matrix::qr_decomp_r() {
    matrix r(col, col);
    matrix q = qr_decomp_q();
    for (int k = 0; k < col; k++) {
        matrix a_k = this->get_col_vec(k);
        for (int i = 0; i <= k; i++) {
            matrix q_i_t = (q.get_col_vec(i)).transpose();
            *r.ref_element(i, k) = (q_i_t * a_k).get_element(0, 0);
        }
    }
    return r;
}
matrix matrix::lu_decomp_l() {
    matrix m(*this);
    matrix l = identity(row);
    for (int k = 0; k < std::min(row, col); k++) {
        int z = k;
        while (!static_cast<bool>(m.arr[z][k])) {
            if (z + 1 == std::min(row, col)) {
                break;
            }
            z++;
        }
        if (k != z) {
            m = m.row_swap(k, z);
        }
        if (m.arr[k][k] == 0) {
            continue;
        }
        // Eliminate entries below the current pivot.
        for (int i = k + 1; i < row; i++) {
            const double multiplier = m.arr[i][k] / m.arr[k][k];
            m = m.row_op(i, -multiplier, k);
            *l.ref_element(i,k) = multiplier;
        }
    }
    neg_zero(l);
    fpg(l);
    return l;
}
matrix matrix::lu_decomp_u() {
    matrix m(*this);
    m.echelon();
    return m;
}

bool matrix::check_upper_tri() {
    // Check the strictly lower-triangular part.
    for (int j = 0; j < col; j++) {
        for (int i = j + 1; i < row; i++) {
            if (get_element(i, j) != 0) {
                return false;
            }
        }
    }
    return true;
}

bool matrix::check_lower_tri() {
    // Check the strictly upper-triangular part.
    for (int i = 0; i < row; i++) {
        for (int j = i + 1; j < col; j++) {
            if (get_element(i, j) != 0) {
                return false;
            }
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
    double c_k = *poly.ref_element(0, 0) = 1;
    // Leverrier-Faddeev iteration builds one coefficient per step.
    for (int i = 1; i <= n; i++) {
        matrix a = (*this) * b_k_1;
        c_k = -a.trace() / i;
        b_k = a + id * c_k;
        *poly.ref_element(0, i) = c_k;
        b_k_1 = b_k;
    }
    return poly;
}

matrix matrix::eigenvalues() {
    if (this->row != this->col) {
        throw std::runtime_error("Eigenvalues are defined only for square matrix.");
    }
    matrix m(*this);
    matrix q = m.qr_decomp_q();
    matrix r = m.qr_decomp_r();
    // Repeated QR steps push the matrix toward upper-triangular form.
    while (!(m.check_upper_tri())) {
        m = r * q;
        q = m.qr_decomp_q();
        r = m.qr_decomp_r();
        fpg(m);
    }
    // The diagonal entries of the converged iterate are the eigenvalue estimates.
    matrix eigen(row, 1);
    for (int i = 0; i < row; i++) {
        *eigen.ref_element(i, 0) = m.get_element(i, i);
    }
    return eigen;
}

double matrix::norm() {
    return sqrt(inner_product(*this, *this));
}

matrix matpow(matrix mat, long long expo) {
    matrix res = identity(mat.row);
    while (expo > 0) {
        if (expo & 1) {
            res = res * mat;
        }
        mat = mat * mat;
        expo >>= 1;
    }
    return res;
}

bool check_ortho(matrix mat) {
    if (mat.transpose().arr == mat.inverse().arr) {
        return true;
    }
    return false;
}

bool check_unitary(matrix mat) {
    if (mat.transpose().arr == mat.adjoint().arr) {
        return true;
    }
    return false;
}

double inner_product(matrix a, matrix b) {
    if ((a.row != b.row) || (a.col != b.col)) {
        throw std::runtime_error("Vectors don't belong same vector space.");
    }
    return ((a.transpose()) * b).trace();
}

double angle(matrix a, matrix b) {
    double l_a = a.norm();
    double l_b = b.norm();
    return acos(inner_product(a, b) / (l_a * l_b));
}

} // namespace linalg
