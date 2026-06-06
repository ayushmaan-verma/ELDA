# ELDA Public Function Index

This file lists the public API declared by the headers under `include/elda/`.

## `include/elda/matrix.hpp`

### Constants

| Name | Description |
| --- | --- |
| `PI` | Pi approximation used by transform helpers and client code. |
| `EPS` | Near-zero threshold used by `fpg()`. |

### `linalg::matrix`

| Signature | Description |
| --- | --- |
| `matrix()` | Construct a `3 x 3` zero matrix. |
| `matrix(int r, int c, double val = 0.0)` | Construct an `r x c` matrix initialized with `val`; negative dimensions throw. |
| `double get_element(int i, int j)` | Return the value at `(i, j)`. |
| `double* ref_element(int i, int j)` | Return a mutable pointer to the value at `(i, j)`. |
| `void input()` | Read every entry from `std::cin`. |
| `void print()` | Print the matrix to `std::cout`. |
| `std::ostream& operator<<(std::ostream& os, const matrix& mat)` | Print the matrix to an output stream. |
| `matrix operator=(matrix m2)` | Assign from a matrix with the same shape. |
| `matrix operator+(matrix m2)` | Add matrices with identical shapes. |
| `matrix operator-(matrix m2)` | Subtract matrices with identical shapes. |
| `matrix operator*(matrix m2)` | Multiply compatible matrices. |
| `matrix operator*(double d)` | Multiply every entry by a scalar. |
| `matrix row_op(int i, double coeff, int j)` | Return a copy with `row[i] += coeff * row[j]`. |
| `matrix col_op(int i, double coeff, int j)` | Return a copy with `col[i] += coeff * col[j]`. |
| `matrix row_swap(int i, int j)` | Return a copy with rows `i` and `j` swapped. |
| `matrix col_swap(int i, int j)` | Return a copy with columns `i` and `j` swapped. |
| `matrix row_multi(int i, double factor)` | Return a copy with row `i` scaled. |
| `matrix col_multi(int i, double factor)` | Return a copy with column `i` scaled. |
| `int echelon()` | Reduce the matrix in place to row-echelon form and return row-swap sign. |
| `int gaussian()` | Reduce the matrix in place with normalized pivots and return row-swap sign. |
| `int gauss_jordan()` | Reduce the matrix in place to reduced row-echelon form. |
| `int canonical()` | Apply row and column elimination toward canonical form. |
| `int rank()` | Return rank computed from Gaussian elimination. |
| `int free_variable()` | Return `row - rank()` for a square coefficient matrix. |
| `double trace()` | Return the sum of the main diagonal of a square matrix. |
| `double diag_product()` | Return the product of entries on the main diagonal. |
| `double det()` | Return the determinant of a square matrix. |
| `double minor(int x, int y)` | Return the minor of entry `(x, y)`. |
| `double cofactor(int x, int y)` | Return the cofactor of entry `(x, y)`. |
| `matrix adjoint()` | Return the classical adjoint, also called the adjugate. |
| `matrix transpose()` | Return the transpose. |
| `matrix inverse()` | Return the inverse of a nonsingular square matrix. |
| `matrix solve()` | Solve an `N x (N + 1)` augmented system and return an `N x 1` solution vector. |
| `matrix orthogonalize()` | Return columns orthogonalized with Gram-Schmidt. |
| `matrix orthonormalize()` | Return columns orthonormalized with Gram-Schmidt. |
| `matrix qr_decomp_q()` | Return the `Q` factor from QR decomposition. |
| `matrix qr_decomp_r()` | Return the `R` factor from QR decomposition. |
| `std::tuple<matrix, matrix, matrix> lu_decomposition()` | Return `{P, L, U}` such that `P * A = L * U`. |
| `matrix get_row(int r)` | Return a matrix containing row `r` and zeros elsewhere. |
| `matrix get_col(int c)` | Return a matrix containing column `c` and zeros elsewhere. |
| `matrix get_row_vec(int r)` | Return row `r` as a `1 x col` matrix. |
| `matrix get_col_vec(int c)` | Return column `c` as a `row x 1` matrix. |
| `void replace_row(int r, matrix rw)` | Replace row `r` with the first row of `rw`. |
| `void replace_col(int c, matrix cn)` | Replace column `c` with the first column of `cn`. |
| `matrix char_poly()` | Return characteristic-polynomial coefficients as a row vector. |
| `matrix eigenvalues()` | Estimate eigenvalues of a square matrix using unshifted QR iteration. |
| `bool check_upper_tri()` | Return true when all strictly lower-triangular entries are zero. |
| `bool check_lower_tri()` | Return true when all strictly upper-triangular entries are zero. |
| `double norm()` | Return the Frobenius norm. |

### Free Matrix Helpers

| Signature | Description |
| --- | --- |
| `bool operator==(matrix m1, matrix m2)` | Return exact entry-wise equality. |
| `bool shape_comp(matrix m1, matrix m2)` | Return true when shapes match. |
| `matrix identity(int n)` | Return an `n x n` identity matrix; negative sizes throw. |
| `void neg_zero(matrix& m)` | Convert exact negative-zero entries to zero. |
| `void fpg(matrix& m)` | Convert entries with absolute value at most `EPS` to zero. |
| `matrix matpow(matrix mat, long long expo)` | Raise a square matrix to a non-negative integer power; throws for negative exponents or non-square inputs. |
| `bool check_ortho(matrix mat)` | Return true when `transpose() == inverse()`. |
| `bool check_unitary(matrix mat)` | Return true when `transpose() == adjoint()`. |
| `double inner_product(matrix a, matrix b)` | Return the Frobenius inner product of same-shaped matrices. |
| `double angle(matrix a, matrix b)` | Return the angle between matrices in radians. |

## `include/elda/transforms.hpp`

| Signature | Description |
| --- | --- |
| `matrix shift(double dx, double dy)` | Return a `3 x 3` homogeneous 2D translation matrix. |
| `matrix scale(double kx, double ky)` | Return a `3 x 3` homogeneous 2D scaling matrix. |
| `matrix rotate(double angle)` | Return a `3 x 3` homogeneous 2D rotation matrix. |
| `matrix shift(double dx, double dy, double dz)` | Return a `4 x 4` homogeneous 3D translation matrix. |
| `matrix scale(double kx, double ky, double kz)` | Return a `4 x 4` homogeneous 3D scaling matrix. |
| `matrix rot_x(double angle)` | Return a `4 x 4` homogeneous rotation matrix around the x-axis. |
| `matrix rot_y(double angle)` | Return a `4 x 4` homogeneous rotation matrix around the y-axis. |
| `matrix rot_z(double angle)` | Return a `4 x 4` homogeneous rotation matrix around the z-axis. |

Angles are measured in radians.

## `include/elda/vector_utils.hpp`

| Signature | Description |
| --- | --- |
| `matrix vec1()` | Return a zero-initialized `1 x 1` column vector. |
| `matrix vec2()` | Return a zero-initialized `2 x 1` column vector. |
| `matrix vec3()` | Return a zero-initialized `3 x 1` column vector. |
| `matrix vec4()` | Return a zero-initialized `4 x 1` column vector. |
| `matrix vec5()` | Return a zero-initialized `5 x 1` column vector. |
| `matrix vec1(double x)` | Return a `1 x 1` column vector containing `x`. |
| `matrix vec2(double x, double y)` | Return a `2 x 1` column vector containing `(x, y)`. |
| `matrix vec3(double x, double y, double z)` | Return a `3 x 1` column vector containing `(x, y, z)`. |
| `matrix vec4(double x, double y, double z, double w)` | Return a `4 x 1` column vector containing `(x, y, z, w)`. |
| `matrix vec5(double x, double y, double z, double w, double u)` | Return a `5 x 1` column vector containing `(x, y, z, w, u)`. |
| `bool check_lin_comb(matrix m1, matrix m2)` | Return true when `m2` is in the span of `m1`. |
| `bool check_lin_comb(matrix m1, matrix m2, matrix m3)` | Return true when `m3` is in the span of `m1` and `m2`. |
| `bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4)` | Return true when `m4` is in the span of `m1`, `m2`, and `m3`. |
| `bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4, matrix m5)` | Return true when `m5` is in the span of `m1`, `m2`, `m3`, and `m4`. |
| `bool check_lin_comb(matrix target, const std::vector<matrix>& vectors)` | Return true when `target` is in the span of the given `vectors`. |
