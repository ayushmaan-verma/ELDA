# ELDA

ELDA is a small C++17 linear algebra library built with CMake. It provides a dense `linalg::matrix` type, matrix operations, a few transformation-matrix helpers, and one vector-style utility helper.

The public headers live under `include/elda/`. The library target is `elda`, and the repository also includes a small demo executable named `main`.

## Features

- Dense matrix storage with public dimensions and element storage
- Matrix addition, subtraction, multiplication, and scalar multiplication
- Elementary row and column operations
- Gaussian, Gauss-Jordan, echelon, and canonical-style reduction helpers
- Determinant, inverse, transpose, adjoint, rank, norm, and characteristic polynomial helpers
- 2D and 3D homogeneous transformation matrix builders

## Repository Layout

```text
lin_alg_lib/
├── CMakeLists.txt
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
    ├── Makefile
    └── elda_booklet.tex
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

This produces:

- `build/libelda.a`
- `build/main`

## Run the Demo

```bash
./build/main
```

The demo constructs a fixed `3 x 3` matrix, applies `gaussian()`, and prints the reduced result.

## Use the Library

Include the full public API:

```cpp
#include <elda/linalg.hpp>
```

Or include only the pieces you need:

```cpp
#include <elda/matrix.hpp>
#include <elda/transforms.hpp>
```

If you consume this repository from another CMake project, link against the `elda` target after adding the source tree:

```cmake
add_subdirectory(path/to/lin_alg_lib)
target_link_libraries(your_target PRIVATE elda)
```

## API Modules

- `matrix.hpp`: the `linalg::matrix` type plus most linear-algebra helpers
- `transforms.hpp`: 2D and 3D homogeneous translation, scaling, and rotation matrices
- `vector_utils.hpp`: `check_similar`, a row-sequence comparison helper
- `linalg.hpp`: umbrella header that includes all public headers

## Behavior Notes

- The library namespace is `linalg`.
- Matrix entries are stored in `matrix::arr` as `std::vector<std::vector<double>>`.
- Many transformation helpers return a new matrix, while reduction helpers such as `echelon()`, `gaussian()`, `gauss_jordan()`, and `canonical()` modify the matrix in place.
- Most dimension mismatches are reported with `std::runtime_error`.
- Angles are interpreted in radians.
- Several low-level helpers assume valid indices and do not perform bounds checking.

## Generate the PDF Booklet

If `latexmk` is installed, build the LaTeX booklet with:

```bash
make -C docs pdf
```

The PDF is written to `docs/build/elda_booklet.pdf`.

## Further Reference

See `DOCUMENTATION.MD` for the detailed API and usage notes.
