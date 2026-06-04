# Contributing to ELDA

Thank you for contributing to ELDA. This guide explains how to make changes that fit the current C++17 library, CMake build, and documentation layout.

## Project Scope

ELDA is a small dense linear algebra library centered on `linalg::matrix`. Contributions should stay focused on:

- matrix operations and decomposition helpers in `include/elda/matrix.hpp` and `src/matrix.cpp`
- homogeneous 2D and 3D transform helpers in `include/elda/transforms.hpp` and `src/transforms.cpp`
- small column-vector helpers in `include/elda/vector_utils.hpp` and `src/vector_utils.cpp`
- project documentation in `README.md`, `DOCUMENTATION.MD`, `FUNCTIONS.MD`, and `docs/`

Avoid broad rewrites unless they are needed for a specific bug fix or feature. Keep public API changes deliberate because users include headers from `include/elda/`.

## Repository Layout

```text
include/elda/        Public headers
src/                 Library implementations
main.cpp             Small demo executable
docs/                Static documentation site
CMakeLists.txt       CMake build for the `elda` library and `main` demo
README.md            Project overview and quick start
DOCUMENTATION.MD     Markdown API documentation
FUNCTIONS.MD         Compact public function index
```

Generated build output belongs in `build/` or another local build directory and should not be committed.

## Development Setup

Requirements:

- CMake 3.14 or newer
- A C++17-compatible compiler such as `g++` or `clang++`
- A native build tool supported by CMake

Configure and build from the repository root:

```bash
cmake -S . -B build
cmake --build build
```

Run the demo executable with a `3 x 3` matrix on standard input:

```bash
printf "3 4 5\n6 7 8\n8 2 3\n" | ./build/main
```

## Coding Guidelines

- Use C++17 and keep code portable across common C++ compilers.
- Put public declarations in `include/elda/` and implementations in the matching file under `src/`.
- Keep all library symbols inside the `linalg` namespace.
- Prefer the existing `matrix` storage model: `row`, `col`, and `arr` backed by `std::vector<std::vector<double>>`.
- Report invalid dimensions with `std::runtime_error`, matching the current API behavior.
- Call `fpg()` when operations can create tiny floating-point residue that should be snapped to zero.
- Preserve the existing split between copy-returning helpers and in-place reduction helpers. Document the behavior clearly when adding a new function.
- Keep comments useful and short. Public header declarations should have concise documentation comments when they are part of the API.
- Do not add unrelated formatting churn to large files.

## API Changes

When adding or changing public API:

- update the relevant header in `include/elda/`
- add or update the implementation in `src/`
- update `FUNCTIONS.MD` with the public signature
- update `DOCUMENTATION.MD` and the static docs in `docs/` when user-facing behavior changes
- consider whether `include/elda/linalg.hpp` should include a new public header

Public behavior notes should mention important constraints such as square-matrix requirements, augmented-matrix shapes, angle units, mutation behavior, and floating-point tolerance.

## Verification

There is currently no formal test suite in the repository. At minimum, run:

```bash
cmake -S . -B build
cmake --build build
printf "3 4 5\n6 7 8\n8 2 3\n" | ./build/main
```

For algorithm changes, also verify small hand-checkable cases. Good examples include:

- identity, zero, diagonal, singular, and non-square matrices
- incompatible dimensions for arithmetic and solving
- determinant, inverse, rank, transpose, and decomposition behavior
- transform helpers with angles in radians
- vector helpers from `vec1()` through `vec5()`

If you add a test framework later, keep tests outside generated build directories and document the command here.

## Documentation Changes

The primary static documentation is `docs/index.html`. If a code change affects documented behavior, update the docs in the same contribution. Keep examples small and executable, and make sure Markdown files and web docs do not disagree about signatures or behavior.

No documentation build step is required. The docs site can be opened directly from `docs/index.html`.

## Commit Checklist

Before submitting a contribution:

- build succeeds with CMake
- demo still runs
- public headers and implementations are in sync
- affected documentation is updated
- generated files from `build/`, `cmake-build-*`, IDE folders, and local artifacts are not committed
- the change is focused on one bug fix, feature, or documentation update


## Running Tests
You can run the automated test suite using CMake and CTest:
```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build