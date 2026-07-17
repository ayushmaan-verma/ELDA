# ELDA

ELDA is a small C++17 linear algebra library built with CMake. It provides a dense `linalg::matrix` type, common matrix operations, decomposition helpers, homogeneous transformation builders, and column-vector utilities.

The public headers live under `include/elda/`. The library target is `elda`, and the repository also includes a small demo executable named `main`.

For the primary documentation set, open `docs/index.html` in a browser. The site is fully static and now organizes the public API by header, with declaration blocks and usage examples for each exported entry, along with behavioral notes and visual 3x3 examples.

## Open Source Contributions

ELDA is open for contributions through GSSoC 2026. Contributors can start by reviewing the [GitHub Issues](https://github.com/ayushmaan-verma/ELDA/issues) for tasks, then follow the guidelines in [`CONTRIBUTING.md`](CONTRIBUTING.md) and the community standards in [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).

For contribution-related questions, contact [ayushmaankumarverma@gmail.com](mailto:ayushmaankumarverma@gmail.com).

## Features

- Dense matrix storage backed by `std::vector<std::vector<double>>`
- Matrix addition, subtraction, multiplication, assignment, and scalar multiplication
- Elementary row and column operations
- Echelon, Gaussian, Gauss-Jordan, and canonical-style reduction helpers
- Determinant, inverse, transpose, adjoint, rank, norm, trace, characteristic polynomial, and linear-system solving helpers
- Gram-Schmidt orthogonalization and orthonormalization
- QR decomposition, LU decomposition, and eigenvalue estimation helpers
- 2D and 3D homogeneous translation, scaling, and rotation matrix builders
- `vec1` through `vec5` helpers for constructing zero or value-filled column vectors
- `check_lin_comb` overloads for span-membership checks on small column-vector sets
- Near-zero floating-point cleanup through `fpg()`

## Repository Layout

```text
lin_alg_lib/
├── CMakeLists.txt
├── main.cpp
├── LICENSE
├── CONTRIBUTING.md
├── CODE_OF_CONDUCT.md
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
├── docs/
│   ├── index.html
│   ├── styles.css
│   └── app.js
├── DOCUMENTATION.MD
└── FUNCTIONS.MD
```

## Requirements

- CMake 3.14 or newer
- A C++17-compatible compiler such as `g++` or `clang++`
- A native build tool supported by CMake

## Build

Configure and compile the library and demo:

```bash
cmake -S . -B build
cmake --build build
```

> **Verification:** Confirm that both `build/libelda.a` and `build/main` are present after the build completes. If either file is missing, rerun the build and check the compiler output for errors.

This produces:

- `build/libelda.a`
- `build/main`

## Testing

Run the automated test suite:

```bash
ctest --test-dir build --output-on-failure
```

## Run the Demo

The demo reads a `3 x 3` matrix from standard input and prints the eigenvalue estimates returned by `matrix::eigenvalues()`.

```bash
printf "3 4 5\n6 7 8\n8 2 3\n" | ./build/main
```

## Use the Library

Include the full public API:

```cpp
#include <elda/linalg.hpp>
```

Or include only the pieces you need:

```cpp
#include <elda/matrix.hpp>
#include <elda/transforms.hpp>
#include <elda/vector_utils.hpp>
```

If you consume this repository from another CMake project, link against the `elda` target after adding the source tree:

```cmake
add_subdirectory(path/to/lin_alg_lib)
target_link_libraries(your_target PRIVATE elda)
```

## API Modules

- `matrix.hpp`: the `linalg::matrix` type plus arithmetic, reductions, decompositions, spectral helpers, and matrix utility functions
- `transforms.hpp`: 2D and 3D homogeneous translation, scaling, and rotation matrices
- `vector_utils.hpp`: `vec1`-`vec5` column-vector constructors and `check_lin_comb` helpers
- `linalg.hpp`: umbrella header that includes all public headers

## Behavior Notes

- The library namespace is `linalg`.
- Matrix entries are stored in `matrix::arr` as `std::vector<std::vector<double>>`.
- Zero-sized matrices are supported, but negative dimensions and negative `identity()` sizes throw `std::runtime_error`.
- Many transformation helpers return a new matrix, while reduction helpers such as `echelon()`, `gaussian()`, `gauss_jordan()`, and `canonical()` modify the matrix in place.
- Most dimension mismatches are reported with `std::runtime_error`.
- Angles are interpreted in radians.
- `trace()`, `det()`, `inverse()`, `char_poly()`, and `eigenvalues()` are defined only for square matrices.
- `solve()` expects an augmented matrix with shape `N x (N + 1)`.
- `EPS` is `1e-6`, and `fpg()` zeros values whose absolute value is at most that threshold. Standard arithmetic operations (`operator+`, `operator-`) and assignment (`operator=`) do not call `fpg()` implicitly, giving the user explicit control over precision cleanup.
- Several low-level helpers assume valid indices and do not perform bounds checking.
- QR and eigenvalue routines use unshifted Gram-Schmidt/QR logic and are intended for simple real-valued workflows.

## Further Reference

- `docs/index.html`: primary static documentation site
- `DOCUMENTATION.MD`: markdown summary aligned with the web docs
- `FUNCTIONS.MD`: compact public function index
- [`CONTRIBUTING.md`](CONTRIBUTING.md): contribution workflow and development guidelines
- [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md): community standards and reporting guidance
- [`LICENSE`](LICENSE): MIT License terms
