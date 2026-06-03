# ELDA Contributor Issue Backlog

This backlog lists issues found during a codebase review. The levels describe the
expected implementation difficulty, not the importance of the bug. Contributors
should add automated tests for every behavior change and update public
documentation when an API contract changes.

## Level 1: Beginner

### 1. Add an automated test target

**Problem:** The repository has no formal test suite or CTest target. Algorithm
changes currently rely on manual checks, which makes regressions difficult to
catch.

**Relevant files:** `CMakeLists.txt`, `CONTRIBUTING.md`, new files under `tests/`

**Acceptance criteria:**

- Add a small C++ test executable and register it with CTest.
- Cover basic construction, arithmetic, transpose, determinant, transforms, and
  vector helpers.
- Document the test command in `CONTRIBUTING.md`.
- `cmake --build build` and `ctest --test-dir build` both succeed.

### 2. Validate matrix dimensions before allocating storage

**Problem:** `matrix(int r, int c)` accepts negative dimensions. These values are
converted to large unsigned sizes by `std::vector`, which can result in an
unexpected allocation exception instead of a clear library error. `identity(n)`
has the same problem for negative `n`.

**Relevant files:** `include/elda/matrix.hpp`, `src/matrix.cpp`

**Acceptance criteria:**

- Reject negative matrix dimensions with a clear `std::runtime_error` or
  `std::invalid_argument`.
- Reject negative identity sizes consistently.
- Add tests for negative rows, negative columns, and negative identity sizes.
- Document whether zero-sized matrices are supported.

### 3. Validate vector shapes in `check_lin_comb`

**Problem:** The `check_lin_comb` overloads assume every argument is a column
vector with the same row count. Mismatched inputs can access `arr[i][0]` outside
the valid range or silently ignore extra rows.

**Relevant files:** `src/vector_utils.cpp`, `include/elda/vector_utils.hpp`

**Acceptance criteria:**

- Verify that every argument has exactly one column.
- Verify that every argument has the same number of rows.
- Throw a clear exception for invalid inputs.
- Add tests for valid spans, mismatched row counts, and non-column matrices.

### 4. Make `matpow` reject unsupported inputs

**Problem:** `matpow` is documented for non-negative powers but silently returns
an identity matrix for a negative exponent. Non-square matrices also produce
confusing behavior, especially for exponent zero.

**Relevant files:** `src/matrix.cpp`, `include/elda/matrix.hpp`,
`DOCUMENTATION.MD`, `FUNCTIONS.MD`, `docs/index.html`

**Acceptance criteria:**

- Reject negative exponents.
- Reject non-square matrices.
- Preserve correct behavior for exponents zero and one.
- Add tests for valid and invalid inputs.

## Level 2: Intermediate

### 5. Fix pivot selection in elimination and rank calculation

**Problem:** `echelon()` and `gaussian()` search for pivots only down to
`min(row, col)` and assume the next pivot is at `(k, k)`. This fails for tall
matrices and for matrices whose leading columns are zero. For example, the
matrix with rows `[0, 1]` and `[0, 2]` is reported as rank 2 even though its
rank is 1.

**Relevant files:** `src/matrix.cpp`

**Acceptance criteria:**

- Track pivot rows and pivot columns independently.
- Search all remaining rows for each candidate pivot column.
- Make `echelon()`, `gaussian()`, `gauss_jordan()`, and `rank()` correct for
  tall, wide, rank-deficient, and leading-zero-column matrices.
- Add regression tests for these matrix shapes.

### 6. Detect singular and inconsistent systems in `solve`

**Problem:** `solve()` assumes every `N x (N + 1)` augmented matrix has one
unique solution. Singular, underdetermined, and inconsistent systems can return
meaningless values instead of reporting that no unique solution exists.

**Relevant files:** `src/matrix.cpp`, `include/elda/matrix.hpp`

**Acceptance criteria:**

- Detect systems with no solution and systems with infinitely many solutions.
- Define and document whether both cases use one exception type or distinct
  error messages.
- Preserve correct solutions for systems that require row swaps.
- Add tests for unique, inconsistent, and underdetermined systems.

### 7. Check for singular pivots before division in `inverse`

**Problem:** `inverse()` divides by the pivot before verifying that it is
non-zero. Singular matrices can create infinities or NaNs during elimination
before the function eventually throws.

**Relevant files:** `src/matrix.cpp`

**Acceptance criteria:**

- Detect a missing pivot before dividing.
- Throw a clear exception for singular matrices without creating non-finite
  intermediate values.
- Use a tolerance consistent with the library's floating-point policy.
- Add tests for zero, duplicate-row, and near-singular matrices.

### 8. Fix Gram-Schmidt column bounds and dependent-column handling

**Problem:** `orthogonalize()` and `orthonormalize()` call
`get_col_vec(i + 1)` on the final loop iteration, which reads past the last
column. They also divide by zero when columns are zero or linearly dependent.

**Relevant files:** `src/matrix.cpp`

**Acceptance criteria:**

- Never request a column index equal to `col`.
- Define behavior for zero and linearly dependent columns.
- Ensure valid outputs contain only finite values.
- Add tests for one-column matrices, rectangular matrices, zero columns, and
  dependent columns.

### 9. Use tolerant comparisons for floating-point matrix checks

**Problem:** `operator==`, `check_ortho()`, and triangular checks use exact
floating-point equality. Matrices produced by trigonometric functions or
decompositions can be mathematically valid but fail these checks due to small
rounding errors.

**Relevant files:** `src/matrix.cpp`, `include/elda/matrix.hpp`

**Acceptance criteria:**

- Define a documented approximate-comparison helper or comparison policy.
- Use it where mathematical checks require tolerance.
- Keep shape mismatches false without accessing invalid storage.
- Add tests using rotation matrices and small floating-point perturbations.

## Level 3: Advanced

### 10. Redesign LU decomposition to handle row permutations

**Problem:** `lu_decomp_l()` swaps rows in its working matrix but does not expose
a permutation matrix or correctly update previously stored multipliers.
`lu_decomp_u()` independently performs elimination. As a result, matrices that
require pivoting do not reliably satisfy `A = L * U`.

**Relevant files:** `src/matrix.cpp`, `include/elda/matrix.hpp`, documentation

**Acceptance criteria:**

- Define an LU API that represents row permutations, such as `P * A = L * U`.
- Keep `L`, `U`, and permutation data consistent through every row swap.
- Document behavior for singular and rectangular matrices.
- Add reconstruction tests for matrices with and without pivoting.

### 11. Add convergence limits and failure handling to `eigenvalues`

**Problem:** `eigenvalues()` loops until the matrix is exactly upper triangular.
Some real matrices have complex eigenvalues, and some QR iterations converge
only asymptotically, so the function can run forever. Exact triangular checks
also make termination sensitive to floating-point residue.

**Relevant files:** `src/matrix.cpp`, `include/elda/matrix.hpp`, documentation

**Acceptance criteria:**

- Add a maximum iteration count and a tolerance-based convergence criterion.
- Report non-convergence clearly instead of looping forever.
- Define behavior for matrices with complex eigenvalues.
- Add tests for diagonal, symmetric, slowly converging, and non-real-spectrum
  cases.

### 12. Make QR decomposition robust for rectangular and rank-deficient inputs

**Problem:** QR decomposition inherits the Gram-Schmidt bounds and zero-division
issues, and its current dimensions and behavior are not clearly defined for
wide or rank-deficient matrices. This also affects `eigenvalues()`.

**Relevant files:** `src/matrix.cpp`, `include/elda/matrix.hpp`, documentation

**Acceptance criteria:**

- Define whether QR returns full or reduced factors for rectangular matrices.
- Ensure factor dimensions are consistent with the documented contract.
- Handle rank-deficient inputs without NaNs or out-of-bounds access.
- Add tests that verify reconstruction and orthogonality within a tolerance.

### 13. Protect matrix invariants from public mutable state

**Problem:** `matrix::row`, `matrix::col`, and `matrix::arr` are public and can
be modified independently. A caller can make the declared shape disagree with
the actual storage, after which most operations can access invalid memory.

**Relevant files:** `include/elda/matrix.hpp`, `src/matrix.cpp`, all public
documentation

**Acceptance criteria:**

- Design a migration path that keeps dimensions and storage consistent.
- Provide safe element access and shape-query APIs.
- Minimize unnecessary source breakage or document the breaking change clearly.
- Add tests showing that public operations cannot observe an inconsistent
  matrix shape.
