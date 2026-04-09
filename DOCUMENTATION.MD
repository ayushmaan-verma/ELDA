# ELDA Documentation

## Overview

ELDA is a C++17 linear algebra library centered around `linalg::matrix`, a dense matrix type backed by `std::vector<std::vector<double>>`. The project builds a static library target named `elda` and a demo executable named `main`.

The public API is defined by:

- `include/elda/linalg.hpp`
- `include/elda/matrix.hpp`
- `include/elda/transforms.hpp`
- `include/elda/vector_utils.hpp`

## Build and Run

Build the library and demo:

```bash
cmake -S . -B build
cmake --build build
```

Run the demo:

```bash
./build/main
```

Build the LaTeX booklet:

```bash
make -C docs pdf
```

## Core Type

### `linalg::matrix`

`matrix` is a mutable dense matrix with three public data members:

- `row`: number of rows
- `col`: number of columns
- `arr`: matrix entries stored row-major as `std::vector<std::vector<double>>`

Constructors:

- `matrix()`: creates a `3 x 3` zero matrix
- `matrix(int r, int c)`: creates an `r x c` zero matrix

Element access:

- `double get_element(int i, int j)`: returns `arr[i][j]`
- `double* ref_element(int i, int j)`: returns a mutable pointer to `arr[i][j]`

Input and output:

- `void input()`: reads all entries from `std::cin`
- `void print()`: writes the matrix to `std::cout`, one row per line

Important notes:

- `get_element()` and `ref_element()` do not perform bounds checking.
- Because `row`, `col`, and `arr` are public, callers can populate matrices directly.

## Constants

Declared in `matrix.hpp`:

- `PI = 3.141593`
- `EPS = 1e-9`

`PI` is used by client code and transformation helpers. `EPS` is declared for floating-point comparisons but is not widely used by the current implementation.

## Arithmetic Operations

Member operators:

- `matrix operator=(matrix m2)`
- `matrix operator+(matrix m2)`
- `matrix operator-(matrix m2)`
- `matrix operator*(matrix m2)`
- `matrix operator*(double d)`

Behavior:

- Assignment requires the left-hand side to already have the same dimensions as the right-hand side.
- Addition and subtraction require identical shapes.
- Matrix multiplication requires `lhs.col == rhs.row`.
- Scalar multiplication returns a new matrix.

All dimension violations throw `std::runtime_error`.

Example:

```cpp
using namespace linalg;

matrix a(2, 2);
matrix b(2, 2);

a.arr = {{1, 2}, {3, 4}};
b.arr = {{5, 6}, {7, 8}};

matrix sum = a + b;
matrix diff = a - b;
matrix prod = a * b;
matrix scaled = a * 2.0;
```

## Shape and Comparison Helpers

Free functions in `matrix.hpp`:

- `bool operator==(matrix m1, matrix m2)`: exact entry-wise equality through `arr`
- `bool shape_comp(matrix m1, matrix m2)`: same `row` and `col`
- `matrix identity(int n)`: returns an `n x n` identity matrix
- `void neg_zero(matrix& m)`: converts exact `-0` values to `0`

These helpers do not apply tolerance-based comparisons.

## Row and Column Operations

These member functions return transformed copies and leave the original matrix unchanged:

- `row_op(i, coeff, j)`: applies `row[i] += coeff * row[j]`
- `col_op(i, coeff, j)`: applies `col[i] += coeff * col[j]`
- `row_swap(i, j)`: swaps rows `i` and `j`
- `col_swap(i, j)`: swaps columns `i` and `j`
- `row_multi(i, factor)`: scales one row
- `col_multi(i, factor)`: scales one column

Example:

```cpp
matrix m(3, 3);
m.arr = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

matrix swapped = m.row_swap(0, 2);
matrix eliminated = m.row_op(1, -4.0, 0);
```

## Reduction and Structural Operations

These member functions modify the matrix in place:

- `int echelon()`
- `int gaussian()`
- `int gauss_jordan()`
- `int canonical()`

Return value:

- Each returns a sign-like integer that tracks row swaps during the forward pass. It starts at `1` and flips when rows are swapped.

Behavior summary:

- `echelon()`: forward elimination without normalizing pivots
- `gaussian()`: forward elimination with each pivot row normalized to `1`
- `gauss_jordan()`: runs `gaussian()` and then clears entries above pivots
- `canonical()`: runs `gaussian()` and then applies column elimination to the right of each pivot

Because these functions mutate the matrix, use a copy if the original data must be preserved.

## Matrix Properties and Derived Matrices

Member functions:

- `int rank()`
- `int free_variable()`
- `double trace()`
- `double diag_product()`
- `double det()`
- `double minor(int x, int y)`
- `double cofactor(int x, int y)`
- `matrix adjoint()`
- `matrix transpose()`
- `matrix inverse()`
- `matrix solve()`
- `matrix get_row(int r)`
- `matrix get_col(int c)`
- `matrix char_poly()`
- `double norm()`

Behavior notes:

- `trace()` and `diag_product()` iterate over `min(row, col)`.
- `det()`, `minor()`, `cofactor()`, `adjoint()`, `inverse()`, and `char_poly()` require a square matrix.
- `free_variable()` requires a square coefficient matrix and returns `row - rank()`.
- `solve()` requires an augmented matrix of shape `N x (N + 1)` and returns an `N x 1` solution vector.
- `get_row(r)` returns a matrix with row `r` preserved and every other row zeroed out.
- `get_col(c)` returns a matrix with column `c` preserved and every other column entry zeroed out.
- `norm()` computes the Frobenius norm.
- `char_poly()` returns a `1 x (n + 1)` row vector of characteristic-polynomial coefficients using the Leverrier-Faddeev recurrence.

Example:

```cpp
matrix m(3, 3);
m.arr = {{3, 4, 5}, {6, 7, 8}, {8, 2, 3}};

double determinant = m.det();
int r = m.rank();
matrix t = m.transpose();
matrix inv = m.inverse();
matrix coeffs = m.char_poly();
```

## Matrix Power and Geometric Helpers

Free functions:

- `matrix matpow(matrix mat, long long expo)`
- `bool check_ortho(matrix mat)`
- `bool check_unitary(matrix mat)`
- `double inner_product(matrix a, matrix b)`
- `double angle(matrix a, matrix b)`

Behavior notes:

- `matpow()` performs repeated squaring and assumes a square matrix with a non-negative exponent.
- `check_ortho()` compares `transpose()` and `inverse()` exactly.
- `check_unitary()` is a library-specific helper that compares `transpose()` and `adjoint()`.
- `inner_product()` computes the Frobenius inner product and requires matching dimensions.
- `angle()` returns the angle in radians from the Frobenius inner product.

Because equality is exact and based on `double` values, `check_ortho()` and `check_unitary()` are sensitive to floating-point roundoff.

## Transformation Helpers

Declared in `transforms.hpp` and implemented as homogeneous transformation matrices.

2D helpers returning `3 x 3` matrices:

- `shift(double dx, double dy)`
- `scale(double kx, double ky)`
- `rotate(double angle)`

3D helpers returning `4 x 4` matrices:

- `shift(double dx, double dy, double dz)`
- `scale(double kx, double ky, double kz)`
- `rot_x(double angle)`
- `rot_y(double angle)`
- `rot_z(double angle)`

Angles are interpreted in radians.

Example:

```cpp
matrix translate2d = shift(2.0, 3.0);
matrix rotate2d = rotate(PI / 4);
matrix scale3d = scale(2.0, 2.0, 2.0);
matrix yaw = rot_y(PI / 6);
```

## Vector Utility Helper

Declared in `vector_utils.hpp`:

- `bool check_similar(matrix mat1, matrix mat2)`

Current behavior:

- removes rows that are exactly all zero
- compares the remaining row sequences for equality

This is not the standard algebraic notion of matrix similarity; it is a library-specific row-pattern helper.

## Error Handling and Limitations

- Dimension mismatches are generally reported with `std::runtime_error`.
- Several functions assume valid row and column indices and do not guard against out-of-range access.
- Exact `double` comparisons are used in many places, so floating-point-heavy workflows may need caller-side tolerance handling.
- `inverse()` is intended for nonsingular square matrices and relies on pivot division during elimination.
- The current CMake setup builds the library and demo locally but does not install or export package configuration files.

## Minimal Example

```cpp
#include <elda/linalg.hpp>

using namespace linalg;

int main() {
    matrix m(2, 2);
    m.arr = {{1, 2}, {3, 4}};

    matrix inv = m.inverse();
    matrix id = m * inv;

    id.print();
}
```
