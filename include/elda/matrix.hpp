#ifndef ELDA_MATRIX_H
#define ELDA_MATRIX_H

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <tuple>

namespace linalg {

/// Pi approximation used by the transform helpers.
constexpr double PI = 3.141593;
/// Threshold used to zero tiny floating-point artifacts.
constexpr double EPS = 1e-6;

/// Dense matrix type backed by a contiguous 1D row-major vector.
class matrix {
  private:
    static size_t validate_dims(int r, int c) {
        if (r < 0 || c < 0) {
            throw std::runtime_error("Matrix dimensions must be non-negative.");
        }
        return static_cast<size_t>(r) * static_cast<size_t>(c);
    }
  public:
    /// Number of rows.
    int row;
    /// Number of columns.
    int col;
    /// Matrix entries flattened in contiguous row-major order.
    std::vector<double> arr;

    /// Constructs a 3x3 zero matrix.
    matrix() : row(3), col(3), arr(9, 0.0) {}

    /// Constructs an r x c zero matrix. Negative dimensions are rejected.
    matrix(int r, int c) : row(r), col(c), arr(validate_dims(r, c), 0.0) {}

    /// Constructs an r x c matrix initialized with val. Negative dimensions are rejected.
    matrix(int r, int c, double val) : row(r), col(c), arr(validate_dims(r, c), val) {}

    // --- Safe Public Element-Access Overloads ---
    double operator()(size_t r, size_t c) const {
        if (r >= static_cast<size_t>(row) || c >= static_cast<size_t>(col)) {
            throw std::out_of_range("Matrix element access out of bounds.");
        }
        return arr[r * col + c];
    }

    double& operator()(size_t r, size_t c) {
        if (r >= static_cast<size_t>(row) || c >= static_cast<size_t>(col)) {
            throw std::out_of_range("Matrix element access out of bounds.");
        }
        return arr[r * col + c];
    }

    /// Returns the number of rows.
    size_t get_rows() const { return static_cast<size_t>(row); }
    /// Returns the number of columns.
    size_t get_cols() const { return static_cast<size_t>(col); }

    /// Returns the value at (i, j) without bounds checking.
    double get_element(int i, int j) { return arr[i * col + j]; }
    /// Returns a mutable pointer to the element at (i, j).
    double* ref_element(int i, int j) { return &arr[i * col + j]; }

    /// Reads matrix entries from `std::cin` in row-major order.
    void input();
    /// Returns a string representation of the matrix.
    std::string to_string() const;
    /// Prints the matrix to `std::cout`, one row per line.
    void print();

    /// Assigns another matrix of the same shape.
    matrix operator=(matrix m2);
    /// Adds two matrices of identical shape.
    matrix operator+(matrix m2);
    /// Subtracts two matrices of identical shape.
    matrix operator-(matrix m2);
    /// Multiplies this matrix by another compatible matrix.
    matrix operator*(matrix m2);
    /// Multiplies every entry by a scalar.
    matrix operator*(double d);
    /// Divides every entry by a scalar.
    matrix operator/(double scalar);

    /// Returns a copy after applying `row[i] += coeff * row[j]`.
    matrix row_op(int i, double coeff, int j);
    /// Returns a copy after applying `col[i] += coeff * col[j]`.
    matrix col_op(int i, double coeff, int j);
    /// Returns a copy with rows `i` and `j` swapped.
    matrix row_swap(int i, int j);
    /// Returns a copy with columns `i` and `j` swapped.
    matrix col_swap(int i, int j);
    /// Returns a copy with row `i` scaled by `factor`.
    matrix row_multi(int i, double factor);
    /// Returns a copy with column `i` scaled by `factor`.
    matrix col_multi(int i, double factor);

    /// Reduces the matrix in place to row-echelon form.
    /// Returns the sign contributed by row swaps.
    int echelon();
    /// Reduces the matrix in place with normalized pivots.
    /// Returns the sign contributed by row swaps.
    int gaussian();
    /// Reduces the matrix in place to reduced row-echelon form.
    /// Returns the sign contributed by row swaps in the forward pass.
    int gauss_jordan();
    /// Applies row and column elimination toward canonical form.
    /// Returns the sign contributed by row swaps in the forward pass.
    int canonical();
    /// Returns the rank computed from Gaussian elimination.
    int rank();
    /// Returns `row - rank()` for square coefficient matrices.
    int free_variable();
    /// Returns the sum of the main diagonal.
    double trace();
    /// Returns the product of the main diagonal.
    double diag_product();
    /// Computes the determinant using echelon reduction.
    double det();
    /// Returns the minor of entry `(x, y)`.
    double minor(int x, int y);
    /// Returns the cofactor of entry `(x, y)`.
    double cofactor(int x, int y);
    /// Returns the classical adjoint (adjugate).
    matrix adjoint();
    /// Returns the transpose.
    matrix transpose();
    /// Returns the inverse using Gauss-Jordan elimination.
    matrix inverse();
    /// Solves an `N x (N + 1)` augmented system and returns the solution column.
    matrix solve();

    /// Orthogonalizes the columns with Gram-Schmidt.
    matrix orthogonalize();
    /// Orthonormalizes the columns with Gram-Schmidt.
    matrix orthonormalize();
    /// Returns the Q factor from QR decomposition.
    matrix qr_decomp_q();
    /// Returns the R factor from QR decomposition.
    matrix qr_decomp_r();

    /// Performs LU decomposition with partial pivoting: P * A = L * U.
    /// Returns a tuple of {P, L, U}. Throws std::runtime_error if matrix is rectangular or singular.
    std::tuple<matrix, matrix, matrix> lu_decomposition() const;

    /// Returns a matrix containing row `r` and zeros elsewhere.
    matrix get_row(int r);
    /// Returns a matrix containing column `c` and zeros elsewhere.
    matrix get_col(int c);
    /// Returns row `r` as a 1 x col vector.
    matrix get_row_vec(int r);
    /// Returns column `c` as a row x 1 vector.
    matrix get_col_vec(int c);
    /// Replaces row `r` with the first row of `rw`.
    void replace_row(int r, matrix rw);
    /// Replaces column `c` with the first column of `cn`.
    void replace_col(int c, matrix cn);

    /// Returns the characteristic polynomial coefficients as a row vector.
    matrix char_poly();
    /// Approximates eigenvalues of a square matrix via unshifted QR iteration.
    /// Returns the diagonal of the converged iterate as a `row x 1` column vector.
    matrix eigenvalues();
    /// Returns true when all strictly lower-triangular entries are exactly zero.
    bool check_upper_tri();
    /// Returns true when all strictly upper-triangular entries are exactly zero.
    bool check_lower_tri();
    /// Returns the Frobenius norm.
    double norm();
};

/// Returns true when both matrices have identical entries.
bool operator==(matrix m1, matrix m2);
/// Returns true when both matrices have the same shape.
bool shape_comp(matrix m1, matrix m2);
/// Returns the `n x n` identity matrix. Negative sizes are rejected.
matrix identity(int n);
/// Normalizes `-0` entries that can appear after elimination.
void neg_zero(matrix& m);
/// Zeros entries whose absolute value is at most `EPS` (floating-point cleanup).
void fpg(matrix& m);
/// Raises a square matrix to a non-negative integer power.
matrix matpow(matrix mat, long long expo);
/// Returns true when `transpose() == inverse()`.
bool check_ortho(matrix mat);
/// Legacy helper that returns true when `transpose() == adjoint()`.
bool check_unitary(matrix mat);
/// Returns the Frobenius inner product of two matrices.
double inner_product(matrix a, matrix b);
/// Returns the angle between two matrices in radians.
double angle(matrix a, matrix b);

std::ostream& operator<<(std::ostream& os, const matrix& m);

}

#endif // ELDA_MATRIX_H