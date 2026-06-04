# ELDA Documentation

## Overview

ELDA is a C++17 linear algebra library centered around `linalg::matrix`, a mutable dense matrix type backed by `std::vector<std::vector<double>>`. The project builds a static library target named `elda` and a demo executable named `main`.

The current primary documentation set is the static web reference at `docs/index.html`. That site is organized by public header and documents every exported entry with its declaration, source mapping, behavior notes, and a usage example.

The public API is defined by:

- `include/elda/linalg.hpp`
- `include/elda/matrix.hpp`
- `include/elda/transforms.hpp`
- `include/elda/vector_utils.hpp`

`linalg.hpp` is the umbrella header and includes the matrix, transform, and vector utility headers.

## Documentation Entry Points

- `docs/index.html`: primary web documentation, organized header-first and function-by-function
- `README.md`: repository overview and build instructions
- `FUNCTIONS.MD`: compact public function index
- `DOCUMENTATION.MD`: repository-local markdown summary aligned with the web docs

Open the web docs directly in a browser:

```bash
xdg-open docs/index.html
```

No documentation build step is required.

## Build and Run

Build the library and demo:

```bash
cmake -S . -B build
cmake --build build
```

Run the demo with nine matrix entries on standard input:

```bash
printf "3 4 5\n6 7 8\n8 2 3\n" | ./build/main
```

The current demo reads a `3 x 3` matrix, computes `matrix::eigenvalues()`, and prints the resulting column vector.

## Header-Oriented API Map

### `include/elda/linalg.hpp`

Umbrella include for the whole public API.

```cpp
#include <elda/linalg.hpp>
using namespace linalg;
```

Use this header when you want matrix algorithms, transforms, and vector helpers in one include.

### `include/elda/matrix.hpp`

Core matrix type plus matrix algorithms.

Primary surface areas:

- constants: `PI`, `EPS`
- public data members: `row`, `col`, `arr`
- construction and I/O: constructors, `input()`, `print()`
- arithmetic: assignment, addition, subtraction, matrix multiplication, scalar multiplication
- elementary operations: row and column combination, swap, and scaling helpers
- reductions: `echelon()`, `gaussian()`, `gauss_jordan()`, `canonical()`
- properties and derived matrices: `rank()`, `trace()`, `det()`, `adjoint()`, `transpose()`, `inverse()`, `solve()`, `norm()`
- structural helpers: `get_row()`, `get_col()`, `get_row_vec()`, `get_col_vec()`, `replace_row()`, `replace_col()`
- decomposition and spectral helpers: `orthogonalize()`, `orthonormalize()`, `qr_decomp_q()`, `qr_decomp_r()`, `lu_decomp_l()`, `lu_decomp_u()`, `char_poly()`, `eigenvalues()`
- free helpers: `identity()`, `fpg()`, `matpow()`, `inner_product()`, `angle()`, and exact-check helpers

Example:

```cpp
#include <elda/matrix.hpp>
using namespace linalg;

matrix a(3, 3);
a.arr = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};

double determinant = a.det();
matrix inverse = a.inverse();
matrix eigen = a.eigenvalues();
```

### `include/elda/transforms.hpp`

Homogeneous transformation matrix builders.

2D helpers return `3 x 3` matrices:

- `shift(double dx, double dy)`
- `scale(double kx, double ky)`
- `rotate(double angle)`

3D helpers return `4 x 4` matrices:

- `shift(double dx, double dy, double dz)`
- `scale(double kx, double ky, double kz)`
- `rot_x(double angle)`
- `rot_y(double angle)`
- `rot_z(double angle)`

Angles are interpreted in radians.

Example:

```cpp
#include <elda/transforms.hpp>
#include <elda/vector_utils.hpp>
using namespace linalg;

matrix translation = shift(2.0, -1.0);
matrix point = vec3(1.0, 2.0, 1.0);
matrix moved = translation * point;
```

### `include/elda/vector_utils.hpp`

Small column-vector factories and span-membership helpers.

Factories:

- `vec1()` through `vec5()`
- `vec1(x)` through `vec5(x, y, z, w, u)`

Span checks:

- `check_lin_comb(matrix m1, matrix m2)`
- `check_lin_comb(matrix m1, matrix m2, matrix m3)`
- `check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4)`
- `check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4, matrix m5)`

Example:

```cpp
#include <elda/vector_utils.hpp>
using namespace linalg;

matrix v1 = vec2(1.0, 2.0);
matrix v2 = vec2(2.0, 4.0);
bool in_span = check_lin_comb(v1, v2);
```

## `linalg::matrix` Summary

### Public Data Members

- `int row`
- `int col`
- `std::vector<std::vector<double>> arr`

The design favors direct manipulation of storage. Most examples populate matrices directly through `arr`.

### Construction

- `matrix()`: constructs a `3 x 3` zero matrix
- `matrix(int r, int c)`: constructs an `r x c` zero matrix

### Element Access

- `get_element(i, j)`: returns `arr[i][j]` without bounds checking
- `ref_element(i, j)`: returns a mutable pointer to `arr[i][j]`

### Mutation Model

ELDA mixes copy-returning helpers and in-place algorithms:

- copy-returning: arithmetic operators, `row_op()`, `row_swap()`, `transpose()`, `inverse()`, `qr_decomp_q()`, `vecN(...)`, and similar helpers
- in-place: `echelon()`, `gaussian()`, `gauss_jordan()`, and `canonical()`

This distinction matters when chaining operations.

## Behavioral Notes

- The library namespace is `linalg`.
- Matrix entries are stored in `matrix::arr` as `std::vector<std::vector<double>>`.
- Most dimension mismatches are reported with `std::runtime_error`.
- Several low-level accessors and replacement helpers assume valid indices and compatible shapes.
- Angles are measured in radians.
- `trace()`, `det()`, `minor()`, `cofactor()`, `adjoint()`, `inverse()`, `char_poly()`, and `eigenvalues()` are defined only for square matrices.
- `solve()` expects an augmented matrix with shape `N x (N + 1)`.
- `EPS` is `1e-6`, and `fpg()` zeros values whose absolute value is at most that threshold.
- Exact `double` comparisons are used in `operator==`, `check_upper_tri()`, `check_lower_tri()`, `check_ortho()`, and `check_unitary()`.
- QR and eigenvalue routines use classical Gram-Schmidt and unshifted QR iteration and are best suited to small real-valued workflows.

## Verified 3x3 Example Outputs

The web documentation includes rendered matrix visuals. The values below match outputs verified from the current codebase.

### Arithmetic Example

Input matrices:

```text
A = [ [ 1,  2, 3],
      [ 0,  1, 4],
      [ 5,  6, 0] ]

B = [ [-2,  1, 0],
      [ 3,  0, 1],
      [ 4, -1, 2] ]
```

Results:

```text
A + B = [ [-1,  3, 3],
          [ 3,  1, 5],
          [ 9,  5, 2] ]

A * B = [ [16, -2, 8],
          [19, -4, 9],
          [ 8,  5, 6] ]
```

### Inverse Example

For

```text
A = [ [1, 2, 3],
      [0, 1, 4],
      [5, 6, 0] ]
```

ELDA returns:

```text
det(A) = 1

A^-1 = [ [-24,  18,  5],
         [ 20, -15, -4],
         [ -5,   4,  1] ]
```

### Linear Solve Example

For the augmented system

```text
[  2   1  -1 |   8 ]
[ -3  -1   2 | -11 ]
[ -2   1   2 |  -3 ]
```

`solve()` returns:

```text
[ 2 ]
[ 3 ]
[-1 ]
```

## Repository Layout

```text
lin_alg_lib/
├── CMakeLists.txt
├── README.md
├── DOCUMENTATION.MD
├── FUNCTIONS.MD
├── main.cpp
├── include/
│   └── elda/
│       ├── linalg.hpp
│       ├── matrix.hpp
│       ├── transforms.hpp
│       └── vector_utils.hpp
├── src/
│   ├── matrix.cpp
│   ├── transforms.cpp
│   └── vector_utils.cpp
└── docs/
    ├── index.html
    ├── styles.css
    └── app.js
```

## Further Reference

- `docs/index.html`: full static web docs
- `FUNCTIONS.MD`: compact API index
- `README.md`: build and repository usage summary

### Eigenvalue Convergence Bounds
The `eigenvalues()` function uses unshifted QR iteration. It enforces a strict iteration limit (1000 iterations) and uses a tolerance-based check (`EPS`) for sub-diagonal convergence. If a matrix has a complex eigenvalue spectrum (e.g., standard rotation matrices), the unshifted real QR algorithm cannot converge to a strictly upper-triangular state. In these asymptotic or non-real scenarios, the function will abort and throw a `std::runtime_error` rather than entering an infinite loop.