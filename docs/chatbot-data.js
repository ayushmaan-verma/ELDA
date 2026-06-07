/**
 * chatbot-data.js
 * ─────────────────────────────────────────────────────────────────────────────
 * ELDA Documentation Chatbot — Knowledge Base
 *
 * Structure of each FAQ entry:
 *   id        — unique string identifier
 *   question  — canonical question text (used for display only)
 *   keywords  — array of lower-case words/phrases used for matching
 *   answer    — plain-text answer shown in the chat bubble
 *   code      — (optional) code snippet shown in a <pre> block
 *   warning   — (optional) important note shown with a ⚠️ badge
 *   note      — (optional) informational note shown with a 💡 badge
 *   related   — array of FAQ ids suggested after this answer
 * ─────────────────────────────────────────────────────────────────────────────
 */

const CHATBOT_DATA = [

  // ── GENERAL ─────────────────────────────────────────────────────────────────

  {
    id: "what_is_elda",
    question: "What is ELDA?",
    keywords: ["what is elda", "about elda", "elda library", "what does elda stand", "introduce elda"],
    answer: "ELDA is a compact C++17 linear algebra library built with CMake. It provides a dense `linalg::matrix` type along with common matrix operations, decomposition helpers, homogeneous transformation builders, and column-vector utilities. It is designed for educational and small real-valued numerical workflows.",
    note: "ELDA is open for contributions through GSSoC 2026.",
    related: ["features", "cpp_version", "is_open_source"]
  },

  {
    id: "what_does_elda_do",
    question: "What does ELDA do?",
    keywords: ["what does elda do", "elda purpose", "elda functionality", "what can elda do"],
    answer: "ELDA lets you work with dense matrices in C++17. It supports arithmetic, row/column operations, elimination routines, decompositions (LU, QR), eigenvalue estimation, characteristic polynomials, 2D/3D homogeneous transforms, and column-vector helpers.",
    related: ["features", "what_is_elda", "which_header"]
  },

  {
    id: "features",
    question: "What are the main features of ELDA?",
    keywords: ["features", "capabilities", "main features", "what features", "what does it support", "functionality"],
    answer: "ELDA's main features:\n• Dense matrix storage backed by std::vector<std::vector<double>>\n• Matrix arithmetic: add, subtract, multiply, scalar multiply\n• Elementary row and column operations\n• Echelon, Gaussian, Gauss-Jordan, and canonical reduction\n• Determinant, inverse, transpose, adjoint, rank, norm, trace\n• Characteristic polynomial and eigenvalue estimation\n• Gram-Schmidt orthogonalization and QR decomposition\n• LU decomposition with partial pivoting\n• Linear system solver\n• 2D and 3D homogeneous translation, scaling, rotation\n• vec1–vec5 column-vector constructors\n• check_lin_comb span-membership checks\n• Near-zero floating-point cleanup via fpg()",
    related: ["what_is_elda", "matrix_ops", "decompositions", "transforms"]
  },

  {
    id: "cpp_version",
    question: "Which C++ version does ELDA support?",
    keywords: ["c++ version", "c++17", "cpp version", "cpp17", "c plus plus", "standard", "compiler"],
    answer: "ELDA targets C++17. You need a C++17-compatible compiler such as g++ or clang++ to build the library.",
    related: ["requirements", "build", "what_is_elda"]
  },

  {
    id: "is_open_source",
    question: "Is ELDA open source?",
    keywords: ["open source", "open-source", "opensource", "free", "free software"],
    answer: "Yes! ELDA is open source. The project is open for contributions, especially through GSSoC 2026. You can find contributor tasks in issues.md.",
    related: ["license", "contribute", "beginner_issues"]
  },

  {
    id: "license",
    question: "What license does ELDA use?",
    keywords: ["license", "mit", "licensing", "terms", "copyright", "usage rights"],
    answer: "ELDA is released under the MIT License. Copyright © 2026 Ayushmaan Kumar Verma. You are free to use, modify, and distribute the software with attribution.",
    related: ["is_open_source", "contribute"]
  },

  // ── INSTALLATION & BUILD ──────────────────────────────────────────────────

  {
    id: "build",
    question: "How do I build ELDA?",
    keywords: ["build", "compile", "cmake", "how to build", "building elda", "how do i build", "setup"],
    answer: "Configure and build the library and demo using CMake. From the repository root, run:",
    code: "cmake -S . -B build\ncmake --build build",
    note: "This produces build/libelda.a (the static library) and build/main (the demo executable).",
    related: ["requirements", "install", "run_demo", "run_tests"]
  },

  {
    id: "requirements",
    question: "What are the requirements?",
    keywords: ["requirements", "prerequisites", "dependencies", "what do i need", "needed", "need to install", "needs"],
    answer: "To build ELDA you need:\n• CMake 3.14 or newer\n• A C++17-compatible compiler (g++ or clang++)\n• A native build tool supported by CMake (e.g. Make, Ninja, MSBuild)",
    related: ["build", "cpp_version"]
  },

  {
    id: "install",
    question: "How do I install ELDA?",
    keywords: ["install", "installation", "how to install", "setup", "add to project", "use in my project"],
    answer: "ELDA is consumed as a CMake subdirectory rather than a traditional install. After building, link the `elda` target in your own CMake project:",
    code: "add_subdirectory(path/to/elda)\ntarget_link_libraries(your_target PRIVATE elda)",
    note: "Include the full public API with: #include <elda/linalg.hpp>",
    related: ["build", "which_header", "requirements"]
  },

  {
    id: "run_demo",
    question: "How do I run the demo?",
    keywords: ["demo", "run demo", "run the demo", "example run", "executable", "main executable"],
    answer: "The demo reads a 3×3 matrix from standard input and prints the eigenvalue estimates. After building, run:",
    code: 'printf "3 4 5\\n6 7 8\\n8 2 3\\n" | ./build/main',
    related: ["build", "eigenvalues", "run_tests"]
  },

  {
    id: "run_tests",
    question: "How do I run tests?",
    keywords: ["run tests", "tests", "testing", "test suite", "unit tests", "how to test"],
    answer: "ELDA uses CTest for automated testing. After building, run the test suite with:",
    code: "ctest --test-dir build --output-on-failure",
    related: ["ctest", "build", "ci"]
  },

  {
    id: "ctest",
    question: "What does ctest do?",
    keywords: ["ctest", "what is ctest", "what does ctest do", "cmake test"],
    answer: "CTest is CMake's built-in test runner. It discovers and runs the test executables registered in CMakeLists.txt. The `--output-on-failure` flag prints the output of any test that fails, making debugging easier.",
    related: ["run_tests", "build"]
  },

  // ── HEADERS ───────────────────────────────────────────────────────────────

  {
    id: "linalg_hpp",
    question: "What is linalg.hpp?",
    keywords: ["linalg.hpp", "linalg hpp", "umbrella header", "linalg header", "main header"],
    answer: "`linalg.hpp` is the umbrella header. Including it brings in the full public API by including matrix.hpp, transforms.hpp, and vector_utils.hpp all at once. It is the recommended entry point for most projects.",
    code: "#include <elda/linalg.hpp>\nusing namespace linalg;",
    related: ["which_header", "matrix_hpp", "transforms_hpp", "vector_utils_hpp"]
  },

  {
    id: "matrix_hpp",
    question: "What is matrix.hpp?",
    keywords: ["matrix.hpp", "matrix hpp", "matrix header"],
    answer: "`matrix.hpp` declares the `linalg::matrix` type along with arithmetic operators, elementary row/column operations, reduction helpers (echelon, gaussian, gauss_jordan, canonical), determinant, inverse, transpose, adjoint, rank, norm, trace, characteristic polynomial, eigenvalue estimation, linear system solver, LU and QR decompositions, and Gram-Schmidt.",
    related: ["linalg_hpp", "which_header", "matrix_ops"]
  },

  {
    id: "transforms_hpp",
    question: "What is transforms.hpp?",
    keywords: ["transforms.hpp", "transforms hpp", "transform header", "transformation header"],
    answer: "`transforms.hpp` provides 2D and 3D homogeneous matrix builders: translation, scaling, and rotation. These return new `linalg::matrix` objects representing the corresponding transformation matrices.",
    related: ["linalg_hpp", "rotation", "scaling", "translation", "radians"]
  },

  {
    id: "vector_utils_hpp",
    question: "What is vector_utils.hpp?",
    keywords: ["vector_utils.hpp", "vector utils hpp", "vector utilities header", "vector header"],
    answer: "`vector_utils.hpp` provides `vec1` through `vec5` constructors for building column vectors of size 1 to 5, and `check_lin_comb` overloads to test whether a vector lies in the span of a given set of column vectors.",
    related: ["linalg_hpp", "vec_helpers", "check_lin_comb"]
  },

  {
    id: "which_header",
    question: "Which header should I include?",
    keywords: ["which header", "what header", "include header", "include file", "what to include"],
    answer: "For most projects, include the umbrella header to get everything:\n\n  #include <elda/linalg.hpp>\n\nIf you only need specific functionality, include selectively:\n• matrix.hpp — matrix type and operations\n• transforms.hpp — 2D/3D transformation matrices\n• vector_utils.hpp — column-vector helpers",
    code: "// Option 1: Full API\n#include <elda/linalg.hpp>\nusing namespace linalg;\n\n// Option 2: Selective\n#include <elda/matrix.hpp>\n#include <elda/transforms.hpp>",
    related: ["linalg_hpp", "matrix_hpp", "transforms_hpp", "vector_utils_hpp"]
  },

  // ── MATRIX OPERATIONS ────────────────────────────────────────────────────

  {
    id: "create_matrix",
    question: "How do I create a matrix?",
    keywords: ["create matrix", "make matrix", "new matrix", "declare matrix", "initialize matrix", "construct matrix"],
    answer: "Construct a `linalg::matrix` with row and column counts, then populate the `arr` member directly. `arr` is a public `std::vector<std::vector<double>>`.",
    code: "linalg::matrix A(3, 3);\nA.arr = {\n    {1, 2, 3},\n    {4, 5, 6},\n    {7, 8, 9}\n};",
    related: ["matrix_storage", "add_matrices", "multiply_matrices"]
  },

  {
    id: "add_matrices",
    question: "How do I add matrices?",
    keywords: ["add matrices", "matrix addition", "plus matrices", "sum matrices", "add two matrices"],
    answer: "Use the `+` operator. It returns a new matrix. Both matrices must have the same dimensions, otherwise a `std::runtime_error` is thrown.",
    code: "linalg::matrix C = A + B;",
    related: ["subtract_matrices", "multiply_matrices", "create_matrix"]
  },

  {
    id: "subtract_matrices",
    question: "How do I subtract matrices?",
    keywords: ["subtract matrices", "matrix subtraction", "minus matrices", "difference matrices"],
    answer: "Use the `-` operator. It returns a new matrix. Both matrices must have the same dimensions.",
    code: "linalg::matrix C = A - B;",
    related: ["add_matrices", "multiply_matrices"]
  },

  {
    id: "multiply_matrices",
    question: "How do I multiply matrices?",
    keywords: ["multiply matrices", "matrix multiplication", "matrix product", "times matrices", "multiply matrix"],
    answer: "Use the `*` operator for matrix–matrix or scalar–matrix multiplication. For matrix multiplication, the number of columns in A must equal the number of rows in B.",
    code: "linalg::matrix C = A * B;   // matrix × matrix\nlinalg::matrix D = 3.0 * A; // scalar × matrix",
    related: ["scalar_mult", "transpose", "add_matrices"]
  },

  {
    id: "scalar_mult",
    question: "How do scalar multiplications work?",
    keywords: ["scalar", "scalar multiply", "scalar multiplication", "multiply by scalar", "scale matrix"],
    answer: "You can multiply a matrix by a `double` scalar using the `*` operator. Scalar multiplication returns a new matrix; it does not modify the original.",
    code: "linalg::matrix B = 2.5 * A; // or A * 2.5",
    related: ["multiply_matrices", "add_matrices"]
  },

  {
    id: "transpose",
    question: "How do I transpose a matrix?",
    keywords: ["transpose", "transposing", "matrix transpose", "t of matrix"],
    answer: "Call the `transpose()` method on a matrix. It returns a new transposed matrix without modifying the original.",
    code: "linalg::matrix At = A.transpose();",
    related: ["adjoint", "inverse", "create_matrix"]
  },

  {
    id: "determinant",
    question: "How do I find the determinant?",
    keywords: ["determinant", "det", "find determinant", "calculate determinant", "matrix determinant"],
    answer: "Call `det()` on a square matrix. It returns a `double`. Throws `std::runtime_error` if the matrix is not square.",
    code: "double d = A.det();",
    warning: "det() is only defined for square matrices.",
    related: ["inverse", "rank", "square_only"]
  },

  {
    id: "inverse",
    question: "How do I find the inverse?",
    keywords: ["inverse", "invert matrix", "matrix inverse", "find inverse", "inverse of a matrix"],
    answer: "Call `inverse()` on a square matrix. It returns a new matrix. Throws if the matrix is not square or is singular.",
    code: "linalg::matrix Ainv = A.inverse();",
    warning: "inverse() requires a square, non-singular matrix.",
    related: ["determinant", "adjoint", "solve_linear", "square_only"]
  },

  {
    id: "trace",
    question: "How do I calculate the trace?",
    keywords: ["trace", "matrix trace", "sum of diagonal", "diagonal sum", "tr(a)"],
    answer: "Call `trace()` on a square matrix. It returns the sum of the main diagonal elements as a `double`.",
    code: "double t = A.trace();",
    warning: "trace() is only defined for square matrices.",
    related: ["determinant", "rank", "square_only"]
  },

  {
    id: "rank",
    question: "How do I calculate the rank?",
    keywords: ["rank", "matrix rank", "find rank", "rank of matrix"],
    answer: "Call `rank()` on a matrix. It returns the rank (number of linearly independent rows/columns) as an `int`.",
    code: "int r = A.rank();",
    related: ["echelon", "gaussian", "determinant"]
  },

  {
    id: "norm",
    question: "How do I calculate the norm?",
    keywords: ["norm", "matrix norm", "frobenius norm", "find norm", "magnitude"],
    answer: "Call `norm()` on a matrix. It returns the Frobenius norm (square root of the sum of squared entries) as a `double`.",
    code: "double n = A.norm();",
    related: ["trace", "rank"]
  },

  {
    id: "solve_linear",
    question: "How do I solve linear equations?",
    keywords: ["solve", "linear equations", "system of equations", "linear system", "solve linear", "ax=b"],
    answer: "Build an augmented matrix with shape N × (N+1) — the coefficient matrix with the RHS column appended — then call `solve()`. It returns the solution vector.",
    code: "// Augmented matrix [A | b]\nlinalg::matrix Ab(3, 4);\nAb.arr = {\n    {2, 1, -1,  8},\n    {-3,-1, 2, -11},\n    {-2, 1, 2, -3}\n};\nlinalg::matrix x = Ab.solve();",
    warning: "solve() expects an augmented matrix of shape N × (N+1).",
    related: ["inverse", "gaussian", "gauss_jordan"]
  },

  // ── REDUCTION METHODS ─────────────────────────────────────────────────────

  {
    id: "echelon",
    question: "What is echelon form?",
    keywords: ["echelon", "row echelon", "echelon form", "ref", "row reduction"],
    answer: "Echelon form (Row Echelon Form) is a matrix form where all zero rows are at the bottom, and the leading entry (pivot) of each non-zero row is to the right of the pivot in the row above. Call `echelon()` to convert in place.",
    code: "A.echelon(); // modifies A in place",
    warning: "echelon() modifies the matrix in place, not returning a new matrix.",
    related: ["gaussian", "gauss_jordan", "canonical", "rank"]
  },

  {
    id: "gaussian",
    question: "What is Gaussian elimination?",
    keywords: ["gaussian", "gaussian elimination", "forward elimination", "gaussian reduction"],
    answer: "`gaussian()` applies Gaussian elimination (forward elimination) to produce row echelon form. It modifies the matrix in place. Used internally for rank computation and solve.",
    code: "A.gaussian(); // modifies A in place",
    warning: "gaussian() is an in-place operation.",
    related: ["echelon", "gauss_jordan", "solve_linear"]
  },

  {
    id: "gauss_jordan",
    question: "What is Gauss-Jordan elimination?",
    keywords: ["gauss jordan", "gauss-jordan", "rref", "reduced row echelon", "gauss jordan elimination"],
    answer: "`gauss_jordan()` applies Gauss-Jordan elimination to produce Reduced Row Echelon Form (RREF). It modifies the matrix in place. Each pivot column is normalized to 1 and all other entries in that column become 0.",
    code: "A.gauss_jordan(); // modifies A in place",
    warning: "gauss_jordan() is an in-place operation.",
    related: ["gaussian", "echelon", "canonical", "solve_linear"]
  },

  {
    id: "canonical",
    question: "What is canonical reduction?",
    keywords: ["canonical", "canonical reduction", "canonical form"],
    answer: "`canonical()` performs a canonical-style reduction on the matrix in place. Like gauss_jordan, it operates in place rather than returning a new matrix.",
    code: "A.canonical(); // modifies A in place",
    warning: "canonical() is an in-place operation.",
    related: ["gauss_jordan", "gaussian", "echelon"]
  },

  // ── DECOMPOSITIONS ────────────────────────────────────────────────────────

  {
    id: "lu_decomp",
    question: "What is LU decomposition?",
    keywords: ["lu decomposition", "lu decomp", "lu", "lu factorization", "plu"],
    answer: "`lu_decomposition()` uses partial pivoting and returns a tuple `{P, L, U}` such that `P * A = L * U`. P is a permutation matrix, L is lower triangular with 1s on the diagonal, and U is upper triangular.",
    code: "auto [P, L, U] = A.lu_decomposition();",
    warning: "lu_decomposition() throws std::runtime_error if the matrix is rectangular or singular.",
    related: ["qr_decomp", "gram_schmidt", "solve_linear"]
  },

  {
    id: "qr_decomp",
    question: "What is QR decomposition?",
    keywords: ["qr decomposition", "qr decomp", "qr", "qr factorization"],
    answer: "`qr_decomposition()` factors the matrix into Q (orthogonal) and R (upper triangular) using classical Gram-Schmidt. Returns a pair `{Q, R}`.",
    code: "auto [Q, R] = A.qr_decomposition();",
    note: "Uses classical (unshifted) Gram-Schmidt — suitable for small, real-valued matrices.",
    related: ["gram_schmidt", "lu_decomp", "eigenvalues"]
  },

  {
    id: "gram_schmidt",
    question: "How does Gram-Schmidt work?",
    keywords: ["gram schmidt", "gram-schmidt", "orthogonalization", "orthonormalize", "orthogonal basis"],
    answer: "ELDA implements classical Gram-Schmidt for QR decomposition and orthogonalization. `gram_schmidt()` orthogonalizes the columns of a matrix; `gram_schmidt_norm()` also normalizes them to unit length.",
    code: "A.gram_schmidt();      // orthogonalize columns\nA.gram_schmidt_norm(); // orthonormalize columns",
    related: ["qr_decomp", "lu_decomp"]
  },

  {
    id: "eigenvalues",
    question: "How are eigenvalues estimated?",
    keywords: ["eigenvalues", "eigenvalue", "eigen values", "spectral", "eigenvalue estimation"],
    answer: "Call `eigenvalues()` on a square matrix. It uses unshifted QR iteration to estimate real eigenvalues and returns them as a column vector (a matrix with 1 column). Convergence is capped at 1000 iterations.",
    code: "linalg::matrix ev = A.eigenvalues();",
    warning: "eigenvalues() throws std::runtime_error if QR iteration fails to converge within 1000 steps. This can happen for matrices with complex eigenvalues.",
    related: ["qr_decomp", "char_poly", "square_only"]
  },

  // ── TRANSFORMATIONS ───────────────────────────────────────────────────────

  {
    id: "translation",
    question: "How do I create translation matrices?",
    keywords: ["translation", "translate", "translation matrix", "homogeneous translation"],
    answer: "Use `translate_2d(tx, ty)` or `translate_3d(tx, ty, tz)` to create homogeneous translation matrices (3×3 and 4×4 respectively).",
    code: "linalg::matrix T2 = linalg::translate_2d(3.0, -1.0);\nlinalg::matrix T3 = linalg::translate_3d(1.0, 2.0, 3.0);",
    related: ["scaling", "rotation", "radians"]
  },

  {
    id: "scaling",
    question: "How do I create scaling matrices?",
    keywords: ["scaling", "scale", "scaling matrix", "homogeneous scaling", "scale matrix"],
    answer: "Use `scale_2d(sx, sy)` or `scale_3d(sx, sy, sz)` to create homogeneous scaling matrices.",
    code: "linalg::matrix S2 = linalg::scale_2d(2.0, 3.0);\nlinalg::matrix S3 = linalg::scale_3d(1.0, 2.0, 0.5);",
    related: ["translation", "rotation", "transforms_hpp"]
  },

  {
    id: "rotation",
    question: "How do I create rotation matrices?",
    keywords: ["rotation", "rotate", "rotation matrix", "homogeneous rotation", "rot_z", "rot_x", "rot_y"],
    answer: "Use `rot_2d(angle)` for a 2D rotation, or `rot_x(angle)`, `rot_y(angle)`, `rot_z(angle)` for 3D rotations around the respective axes. All angles are in radians.",
    code: "#include <cmath>\nlinalg::matrix R2 = linalg::rot_2d(M_PI / 4); // 45°\nlinalg::matrix Rx = linalg::rot_x(M_PI / 2); // 90° around X",
    warning: "All transformation angles are interpreted in radians, not degrees.",
    related: ["translation", "scaling", "radians"]
  },

  {
    id: "radians",
    question: "Are angles in radians?",
    keywords: ["radians", "degrees", "angle units", "angle in radians", "angles"],
    answer: "Yes. All rotation and transformation angles in ELDA are interpreted in radians. Convert degrees to radians using: radians = degrees × π / 180.",
    code: "// Convert 45 degrees to radians\ndouble angle = 45.0 * M_PI / 180.0;\nlinalg::matrix R = linalg::rot_2d(angle);",
    related: ["rotation", "scaling", "translation"]
  },

  // ── VECTORS ───────────────────────────────────────────────────────────────

  {
    id: "vec_helpers",
    question: "What are vec1 to vec5?",
    keywords: ["vec1", "vec2", "vec3", "vec4", "vec5", "vector helpers", "column vector helpers", "vec helper"],
    answer: "`vec1()` through `vec5()` are convenience constructors that create column vectors (n×1 matrices) of size 1 to 5. Pass values directly as arguments, or pass a single value to fill all entries with that value.",
    code: "linalg::matrix v = linalg::vec3(1.0, 2.0, 3.0);\n// Creates a 3×1 column vector [1, 2, 3]^T",
    related: ["create_column_vec", "check_lin_comb", "vector_utils_hpp"]
  },

  {
    id: "create_column_vec",
    question: "How do I create column vectors?",
    keywords: ["column vector", "create column vector", "column vec", "make vector", "create vector"],
    answer: "Use the `vec1`–`vec5` helpers for small vectors, or manually create a matrix with 1 column (`linalg::matrix(n, 1)`) and populate `arr`.",
    code: "// Using helper\nlinalg::matrix v = linalg::vec3(1.0, 0.0, -1.0);\n\n// Manual\nlinalg::matrix w(4, 1);\nw.arr = {{2.0}, {3.0}, {-1.0}, {0.0}};",
    related: ["vec_helpers", "check_lin_comb"]
  },

  {
    id: "check_lin_comb",
    question: "What is check_lin_comb?",
    keywords: ["check_lin_comb", "linear combination", "span check", "linear combo", "lin comb"],
    answer: "`check_lin_comb()` tests whether a target column vector can be expressed as a linear combination of a given set of column vectors (i.e., whether it lies in their span). Returns `true` if it can.",
    code: "linalg::matrix a = linalg::vec2(1.0, 0.0);\nlinalg::matrix b = linalg::vec2(0.0, 1.0);\nlinalg::matrix v = linalg::vec2(3.0, 4.0);\nbool result = linalg::check_lin_comb(v, {a, b}); // true",
    warning: "Inputs must be strictly column vectors (cols == 1) with identical row counts. Mismatches throw std::runtime_error.",
    related: ["vec_helpers", "span_check"]
  },

  {
    id: "span_check",
    question: "How do span checks work?",
    keywords: ["span", "span check", "in span", "linearly independent", "basis check"],
    answer: "`check_lin_comb()` performs a span check by forming an augmented matrix from the given vectors and the target, then applying Gauss-Jordan to test for consistency. It throws if vectors are not column vectors or have inconsistent dimensions.",
    related: ["check_lin_comb", "vec_helpers", "gauss_jordan"]
  },

  // ── BEHAVIOR NOTES ────────────────────────────────────────────────────────

  {
    id: "namespace",
    question: "What namespace does ELDA use?",
    keywords: ["namespace", "linalg namespace", "using namespace", "name space"],
    answer: "All ELDA symbols live in the `linalg` namespace. Use `using namespace linalg;` for convenience, or qualify each type/function as `linalg::matrix`, `linalg::translate_2d`, etc.",
    code: "#include <elda/linalg.hpp>\nusing namespace linalg; // optional",
    related: ["what_is_elda", "which_header"]
  },

  {
    id: "matrix_storage",
    question: "How are matrices stored?",
    keywords: ["storage", "matrix storage", "how stored", "internal representation", "arr", "matrix arr"],
    answer: "Matrix entries are stored in the public member `matrix::arr`, which is a `std::vector<std::vector<double>>`. The outer vector indexes rows; each inner vector holds the columns of that row. `matrix::row` and `matrix::col` store the dimensions.",
    code: "linalg::matrix A(2, 3); // 2 rows, 3 cols\nA.arr[0][1] = 5.0;       // row 0, col 1",
    related: ["create_matrix", "namespace"]
  },

  {
    id: "eps",
    question: "What is EPS?",
    keywords: ["eps", "epsilon", "eps constant", "tolerance", "floating point tolerance"],
    answer: "`EPS` is a constant equal to `1e-6`. It serves as the near-zero threshold used by `fpg()` to clean up floating-point residue after arithmetic operations.",
    related: ["fpg", "matrix_storage"]
  },

  {
    id: "fpg",
    question: "What does fpg() do?",
    keywords: ["fpg", "fpg()", "floating point guard", "zero cleanup", "snap to zero", "near zero"],
    answer: "`fpg()` zeros out any matrix entry whose absolute value is at or below `EPS` (1e-6). Call it after arithmetic or elimination operations to clean up tiny floating-point residues that should be exactly zero.",
    code: "A.fpg(); // snaps near-zero values to 0.0",
    related: ["eps", "matrix_storage"]
  },

  {
    id: "square_only",
    question: "Which operations require square matrices?",
    keywords: ["square matrix", "square only", "non-square", "requires square", "square matrices"],
    answer: "The following operations require square matrices and will throw `std::runtime_error` otherwise:\n• trace()\n• det()\n• minor()\n• cofactor()\n• adjoint()\n• inverse()\n• char_poly()\n• eigenvalues()\n• lu_decomposition() (also rejects rectangular)",
    warning: "Calling these on non-square matrices throws std::runtime_error.",
    related: ["determinant", "inverse", "trace", "eigenvalues", "exceptions"]
  },

  {
    id: "exceptions",
    question: "What exceptions can ELDA throw?",
    keywords: ["exception", "exceptions", "throw", "error", "runtime error", "std::runtime_error", "what can throw"],
    answer: "ELDA reports errors using `std::runtime_error`. Common causes:\n• Negative dimensions or negative identity() size\n• Dimension mismatch in arithmetic\n• Square-only operations called on non-square matrices\n• Singular matrix passed to inverse() or lu_decomposition()\n• solve() called with wrong augmented matrix shape (not N×N+1)\n• eigenvalues() failing to converge within 1000 iterations\n• check_lin_comb() with non-column vectors or mismatched sizes",
    related: ["square_only", "solve_linear", "lu_decomp", "eigenvalues"]
  },

  // ── DOCUMENTATION ─────────────────────────────────────────────────────────

  {
    id: "examples",
    question: "Where can I find examples?",
    keywords: ["examples", "find examples", "usage examples", "code examples", "sample code", "demo"],
    answer: "Examples are available in multiple places:\n• docs/index.html — visual 3×3 worked examples (Arithmetic, Inverse, Linear Solve, QR, 2D Transforms)\n• DOCUMENTATION.MD — markdown API documentation with snippets\n• main.cpp — the demo executable source",
    related: ["documentation_md", "functions_md", "run_demo"]
  },

  {
    id: "where_docs",
    question: "Where is the main documentation?",
    keywords: ["documentation", "where is docs", "main docs", "primary docs", "find documentation", "docs site"],
    answer: "The primary documentation is the static site at `docs/index.html`. Open it directly in any browser — no build step required. It covers the full public API organized by header, with examples and behavior notes.",
    related: ["documentation_md", "functions_md", "examples"]
  },

  {
    id: "documentation_md",
    question: "What is DOCUMENTATION.MD?",
    keywords: ["documentation.md", "documentation md", "docs md", "markdown docs"],
    answer: "`DOCUMENTATION.MD` is a markdown summary of the ELDA API aligned with the web documentation site. It provides a human-readable reference without needing a browser.",
    related: ["functions_md", "where_docs"]
  },

  {
    id: "functions_md",
    question: "What is FUNCTIONS.MD?",
    keywords: ["functions.md", "functions md", "function index", "api index", "function list"],
    answer: "`FUNCTIONS.MD` is a compact public function index listing every exported function signature in ELDA. It is a quick reference for looking up exact function names and signatures.",
    related: ["documentation_md", "where_docs"]
  },

  // ── CONTRIBUTING ──────────────────────────────────────────────────────────

  {
    id: "contribute",
    question: "How can I contribute?",
    keywords: ["contribute", "contributing", "how to contribute", "pull request", "pr", "open source contribution", "submit"],
    answer: "To contribute to ELDA:\n1. Review issues.md for tasks grouped by difficulty\n2. Follow the guidelines in CONTRIBUTING.md\n3. Keep all library symbols in the linalg namespace\n4. Use C++17 and match the existing coding style\n5. Build and pass tests before submitting\n6. Update relevant documentation alongside code changes",
    note: "ELDA is participating in GSSoC 2026 — great time to contribute!",
    related: ["beginner_issues", "contributing_md", "code_of_conduct", "contact"]
  },

  {
    id: "beginner_issues",
    question: "Where are beginner issues?",
    keywords: ["beginner issues", "beginner", "easy issues", "good first issue", "starter issues", "newbie", "first issue", "beginner tasks"],
    answer: "Beginner-friendly tasks are listed in `issues.md`, grouped by difficulty level. Look for tasks tagged as easy or beginner. The file is in the repository root.",
    related: ["contribute", "contributing_md"]
  },

  {
    id: "contributing_md",
    question: "What is CONTRIBUTING.md?",
    keywords: ["contributing.md", "contributing md", "contribution guide", "contribution guidelines"],
    answer: "`CONTRIBUTING.md` documents the full contribution workflow: development setup, coding guidelines, API change process, verification steps (build, lint with cppcheck, format with clang-format), and a commit checklist.",
    related: ["contribute", "code_of_conduct", "run_tests"]
  },

  {
    id: "code_of_conduct",
    question: "What is CODE_OF_CONDUCT.md?",
    keywords: ["code of conduct", "code_of_conduct.md", "conduct", "community standards", "behavior policy"],
    answer: "`CODE_OF_CONDUCT.md` outlines the community standards for the ELDA project. It covers expected behavior, reporting guidelines, and consequences for violations. All contributors are expected to follow it.",
    related: ["contributing_md", "contribute"]
  },

  {
    id: "contact",
    question: "How can I contact the maintainer?",
    keywords: ["contact", "maintainer", "email", "reach out", "owner", "contact maintainer", "report issue"],
    answer: "For contribution-related questions, contact the project owner Ayushmaan Kumar Verma at:\n📧 ayushmaankumarverma@gmail.com",
    related: ["contribute", "contributing_md", "code_of_conduct"]
  },

  {
    id: "greetings",
    question: "Greetings",
    keywords: [
      "hi", "hello", "hey", "hii", "hlo", 
      "good morning", "good afternoon", "good evening", 
      "namaste", "greetings", "yo", "sup"
    ],
    answer: "Hello! 👋 I'm ELDA Assistant. I can help you explore the ELDA documentation, explain matrix operations, guide you through installation, and answer questions about the library. Try asking: 'How do I build ELDA?' or 'What is matrix.hpp?'",
    related: ["what_is_elda", "build", "matrix_hpp"]
  },

  {
    id: "thanks",
    question: "Thank you",
    keywords: ["thanks", "thank you"],
    answer: "You're very welcome! 😊 Let me know if you have any other questions about ELDA or if there's anything else I can help you with.",
    related: ["what_is_elda", "contribute"]
  },

  {
    id: "goodbye",
    question: "Goodbye",
    keywords: ["bye", "goodbye", "see you"],
    answer: "Goodbye! 👋 Happy coding with ELDA! Don't hesitate to reach out if you need more help later.",
    related: ["what_is_elda"]
  },

  {
    id: "help_info",
    question: "How can you help me?",
    keywords: ["who are you", "what can you do", "help", "can you help me", "how can you help me"],
    answer: "I am **ELDA Assistant**, your documentation helper. I can:\n• Guide you through **installation** and building\n• Explain **headers** like `matrix.hpp`, `transforms.hpp`, etc.\n• Provide examples of **matrix operations** (det, inverse, transpose)\n• Explain **LU/QR decompositions** and eigenvalue estimation\n• Guide you on how to **contribute** for GSSoC 2026\n\nTry asking: 'How do I build ELDA?' or 'What is matrix.hpp?'",
    related: ["what_is_elda", "build", "contribute"]
  }

];

/**
 * Bonus: Suggested "getting started" questions shown in the welcome message.
 * These are displayed as quick-reply chips to guide new users.
 */
const STARTER_QUESTIONS = [
  "What is ELDA?",
  "How do I build ELDA?",
  "How do I create a matrix?",
  "What are the main features?",
  "How can I contribute?"
];
