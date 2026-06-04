#ifndef ELDA_MATRIX_H
#define ELDA_MATRIX_H

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

namespace linalg {

/// Pi approximation used by the transform helpers.
constexpr double PI = 3.141593;
/// Threshold used to zero tiny floating-point artifacts.
constexpr double EPS = 1e-6;

/// Dense matrix type backed by a contiguous 1D row-major vector.
class matrix {
  public:
    /// Number of rows.
    int row;
    /// Number of columns.
    int col;
    /// Matrix entries flattened in contiguous row-major order.
    std::vector<double> arr;

    /// Constructs a 3x3 zero matrix using modern STL block initialization.
    matrix() : row(3), col(3), arr(9, 0.0) {}

    /// Constructs an r x c zero matrix using modern STL block initialization.
    matrix(int r, int c) : row(r), col(c), arr(r * c, 0.0) {}

    /// Returns true if the matrix is square (rows == columns).
    bool is_square() const { return row == col; }

    inline int get_rows() const { return row; }
    inline int get_cols() const { return col; }

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

    /// Returns the value at (i, j) without bounds checking.
    double get_element(int i, int j) { return arr[i * col + j]; }
    /// Returns a mutable pointer to the element at (i, j).
    double* ref_element(int i, int j) { return &arr[i * col + j]; }

    /// Returns a string representation of the matrix.
    std::string to_string() const;

    /// Reads matrix entries from std::cin in row-major order.
    void input();
    /// Prints the matrix to std::cout, one row per line.
    void print();

    /// Assigns another matrix of the same shape.
    matrix operator=(const matrix& m2);
    /// Adds two matrices of identical shape.
    matrix operator+(const matrix& m2);
    /// Subtracts two matrices of identical shape.
    matrix operator-(const matrix& m2);
    /// Multiplies this matrix by another compatible matrix.
    matrix operator*(const matrix& m2);
    /// Multiplies every entry by a scalar.
    matrix operator*(double d);

    matrix row_op(int i, double coeff, int j);
    matrix col_op(int i, double coeff, int j);
    matrix row_swap(int i, int j);
    matrix col_swap(int i, int j);
    matrix row_multi(int i, double factor);
    matrix col_multi(int i, double factor);

    int echelon();
    int gaussian();
    int gauss_jordan();
    int canonical();
    int rank();
    int free_variable();
    double trace();
    double diag_product();
    double det();
    double minor(int x, int y);
    double cofactor(int x, int y);
    matrix adjoint();
    matrix transpose();
    matrix inverse();
    matrix solve();

    matrix orthogonalize();
    matrix orthonormalize();
    matrix qr_decomp_q();
    matrix qr_decomp_r();
    std::tuple<matrix, matrix, matrix> lu_decomposition() const;

    matrix get_row(int r);
    matrix get_col(int c);
    matrix get_row_vec(int r);
    matrix get_col_vec(int c);
    void replace_row(int r, const matrix& rw);
    void replace_col(int c, const matrix& cn);

    matrix char_poly();
    matrix eigenvalues();
    bool check_upper_tri();
    bool check_lower_tri();
    double norm();
};

bool operator==(const matrix& m1, const matrix& m2);
bool shape_comp(const matrix& m1, const matrix& m2);
matrix identity(int n);
void neg_zero(matrix& m);
void fpg(matrix& m);
matrix matpow(matrix mat, long long expo);
bool check_ortho(const matrix& mat);
bool check_unitary(const matrix& mat);
double inner_product(matrix a, matrix b);
double angle(matrix a, matrix b);

}

#endif // ELDA_MATRIX_H