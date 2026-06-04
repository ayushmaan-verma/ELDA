function createEntry(config) {
    return {
        kind: "function",
        implementedIn: "see source",
        notes: [],
        ...config
    };
}

const referenceData = [
    {
        id: "linalg",
        mountId: "api-linalg",
        title: "Umbrella include",
        headerPath: "include/elda/linalg.hpp",
        synopsis: `#include <elda/linalg.hpp>`,
        intro:
            "This header is the one-stop entry point for consumers who want the entire public ELDA surface in one include.",
        chips: ["1 export surface", "umbrella header"],
        groups: [
            {
                title: "Exported API",
                copy: "The umbrella header re-exports the three concrete public headers.",
                entries: [
                    createEntry({
                        id: "linalg-include",
                        title: "Umbrella include",
                        signature: `#include <elda/linalg.hpp>`,
                        kind: "header",
                        declaredIn: "include/elda/linalg.hpp",
                        implementedIn: "header-only include wrapper",
                        description:
                            "Includes matrix, transform, and vector utility declarations in a single directive.",
                        example: `#include <elda/linalg.hpp>
using namespace linalg;

matrix a(3, 3);
matrix t = shift(1.0, 2.0);
matrix v = vec3(1.0, 0.0, 1.0);`,
                        notes: [
                            "Equivalent to including matrix.hpp, transforms.hpp, and vector_utils.hpp manually.",
                            "Recommended for examples and small applications."
                        ]
                    })
                ]
            }
        ]
    },
    {
        id: "matrix",
        mountId: "api-matrix",
        title: "Core matrix API",
        headerPath: "include/elda/matrix.hpp",
        synopsis: `#include <elda/matrix.hpp>`,
        intro:
            "This header declares the linalg::matrix type, its algorithms, and the free helpers that operate on matrices.",
        chips: ["dense matrices", "algorithms", "free helpers"],
        groups: [
            {
                title: "Constants",
                copy: "Numerical constants used by the matrix and transform layers.",
                entries: [
                    createEntry({
                        id: "matrix-pi",
                        title: "PI",
                        signature: `constexpr double PI = 3.141593;`,
                        kind: "constant",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "header constant",
                        description: "Pi approximation exposed for client code and transform helpers.",
                        example: `double quarter_turn = PI / 2.0;
matrix r = rotate(quarter_turn);`,
                        notes: ["Used throughout examples that require angle values in radians."]
                    }),
                    createEntry({
                        id: "matrix-eps",
                        title: "EPS",
                        signature: `constexpr double EPS = 1e-6;`,
                        kind: "constant",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "header constant",
                        description: "Near-zero threshold used by fpg() to suppress floating-point residue.",
                        example: `if (std::abs(value) <= EPS) {
    value = 0.0;
}`,
                        notes: ["The cleanup helper fpg() applies this threshold element-wise."]
                    })
                ]
            },
            {
                title: "Public data members",
                copy: "The matrix class exposes its shape metadata and storage directly.",
                entries: [
                    createEntry({
                        id: "matrix-row-member",
                        title: "row",
                        signature: `int row;`,
                        kind: "data member",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "stored in each matrix instance",
                        description: "Number of rows in the matrix.",
                        example: `matrix m(3, 4);
int rows = m.row;`,
                        notes: ["This value is initialized by the constructors and used throughout the algorithms."]
                    }),
                    createEntry({
                        id: "matrix-col-member",
                        title: "col",
                        signature: `int col;`,
                        kind: "data member",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "stored in each matrix instance",
                        description: "Number of columns in the matrix.",
                        example: `matrix m(3, 4);
int cols = m.col;`,
                        notes: ["This value is public and can be read directly by callers."]
                    }),
                    createEntry({
                        id: "matrix-arr-member",
                        title: "arr",
                        signature: `std::vector<std::vector<double>> arr;`,
                        kind: "data member",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "stored in each matrix instance",
                        description: "Row-major storage for the matrix entries.",
                        example: `matrix m(2, 2);
m.arr = {{1.0, 2.0}, {3.0, 4.0}};`,
                        notes: ["Most examples in the documentation populate matrices directly through arr."]
                    })
                ]
            },
            {
                title: "Construction, element access, and I/O",
                copy: "The class is mutable and exposes its storage publicly.",
                entries: [
                    createEntry({
                        id: "matrix-default-ctor",
                        title: "Default constructor",
                        signature: `matrix()`,
                        kind: "constructor",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "inline in header",
                        description: "Constructs a 3 x 3 zero matrix.",
                        example: `matrix m;
m.print();`,
                        notes: ["The default shape is fixed at 3 x 3.", "Storage is initialized to zero."]
                    }),
                    createEntry({
                        id: "matrix-sized-ctor",
                        title: "Sized constructor",
                        signature: `matrix(int r, int c, double val = 0.0)`,
                        kind: "constructor",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "inline in header",
                        description: "Constructs an r x c matrix initialized with val.",
                        example: `matrix m(2, 2, 1.0);`,
                        notes: ["Storage is allocated to r by c and every entry is initialized to val.", "val defaults to 0.0.", "Zero-sized matrices are supported.", "Negative dimensions throw before allocation."]
                    }),
                    createEntry({
                        id: "matrix-get-element",
                        title: "get_element",
                        signature: `double get_element(int i, int j)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "inline in header",
                        description: "Returns arr[i][j] without bounds checking.",
                        example: `matrix m(2, 2);
m.arr = {{1.0, 2.0}, {3.0, 4.0}};
double value = m.get_element(1, 0);`,
                        notes: ["Callers are responsible for providing valid indices."]
                    }),
                    createEntry({
                        id: "matrix-ref-element",
                        title: "ref_element",
                        signature: `double* ref_element(int i, int j)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "inline in header",
                        description: "Returns a mutable pointer to arr[i][j].",
                        example: `matrix m(2, 2);
*m.ref_element(0, 1) = 7.5;`,
                        notes: ["Useful in the transform helpers and direct manual updates."]
                    }),
                    createEntry({
                        id: "matrix-input",
                        title: "input",
                        signature: `void input()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Reads matrix elements from std::cin in row-major order.",
                        example: `matrix m(2, 2);
m.input();   // expected input: four numbers`,
                        notes: ["The function performs no prompt or format validation."]
                    }),
                    createEntry({
                        id: "matrix-print",
                        title: "print",
                        signature: `void print()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Prints the matrix row by row to std::cout.",
                        example: `matrix m(2, 2);
m.arr = {{1.0, 2.0}, {3.0, 4.0}};
m.print();`,
                        notes: ["Entries are printed with spaces between columns."]
                    })
                ]
            },
            {
                title: "Operators and core helpers",
                copy: "Arithmetic is dense and eager, with near-zero cleanup on most return values.",
                entries: [
                    createEntry({
                        id: "matrix-assign",
                        title: "operator=",
                        signature: `matrix operator=(matrix m2)`,
                        kind: "operator",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Copies values from another matrix when both operands already have the same shape.",
                        example: `matrix lhs(2, 2);
matrix rhs(2, 2);
rhs.arr = {{1.0, 2.0}, {3.0, 4.0}};

lhs = rhs;`,
                        notes: [
                            "This custom operator does not resize the left-hand side.",
                            "A shape mismatch throws std::runtime_error."
                        ]
                    }),
                    createEntry({
                        id: "matrix-add",
                        title: "operator+",
                        signature: `matrix operator+(matrix m2)`,
                        kind: "operator",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Adds two matrices of identical shape and returns a new matrix.",
                        example: `matrix a(2, 2);
matrix b(2, 2);
a.arr = {{1.0, 2.0}, {3.0, 4.0}};
b.arr = {{5.0, 6.0}, {7.0, 8.0}};

matrix sum = a + b;`,
                        notes: ["The result is normalized with fpg() before return."]
                    }),
                    createEntry({
                        id: "matrix-subtract",
                        title: "operator-",
                        signature: `matrix operator-(matrix m2)`,
                        kind: "operator",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Subtracts two matrices of identical shape and returns a new matrix.",
                        example: `matrix a(2, 2);
matrix b(2, 2);
a.arr = {{5.0, 6.0}, {7.0, 8.0}};
b.arr = {{1.0, 2.0}, {3.0, 4.0}};

matrix diff = a - b;`,
                        notes: ["A shape mismatch throws std::runtime_error."]
                    }),
                    createEntry({
                        id: "matrix-multiply",
                        title: "operator* (matrix)",
                        signature: `matrix operator*(matrix m2)`,
                        kind: "operator",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Performs dense matrix multiplication when the dimensions are compatible.",
                        example: `matrix a(2, 3);
matrix b(3, 1);
a.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
b.arr = {{1.0}, {0.0}, {1.0}};

matrix product = a * b;`,
                        notes: ["Requires lhs.col == rhs.row."]
                    }),
                    createEntry({
                        id: "matrix-scale-operator",
                        title: "operator* (scalar)",
                        signature: `matrix operator*(double d)`,
                        kind: "operator",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Scales every element by a scalar and returns a new matrix.",
                        example: `matrix a(2, 2);
a.arr = {{1.0, -2.0}, {3.0, 4.0}};

matrix scaled = a * 0.5;`,
                        notes: ["The matrix itself is not modified."]
                    }),
                    createEntry({
                        id: "matrix-equality",
                        title: "operator==",
                        signature: `bool operator==(matrix m1, matrix m2)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Checks exact entry-wise equality through the underlying arr storage.",
                        example: `matrix a(2, 2);
matrix b(2, 2);
a.arr = {{1.0, 2.0}, {3.0, 4.0}};
b.arr = {{1.0, 2.0}, {3.0, 4.0}};

bool same = (a == b);`,
                        notes: ["This is not a tolerance-based numerical comparison."]
                    }),
                    createEntry({
                        id: "matrix-shape-comp",
                        title: "shape_comp",
                        signature: `bool shape_comp(matrix m1, matrix m2)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns true when two matrices have the same row and column counts.",
                        example: `matrix a(3, 2);
matrix b(3, 2);
bool same_shape = shape_comp(a, b);`,
                        notes: ["Only shape is checked, not the stored values."]
                    }),
                    createEntry({
                        id: "matrix-identity",
                        title: "identity",
                        signature: `matrix identity(int n)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Creates an n x n identity matrix.",
                        example: `matrix i3 = identity(3);`,
                        notes: ["The diagonal is filled with 1 and the remaining entries stay 0.", "Negative sizes throw before allocation."]
                    }),
                    createEntry({
                        id: "matrix-neg-zero",
                        title: "neg_zero",
                        signature: `void neg_zero(matrix& m)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Rewrites exact negative-zero entries to ordinary zero.",
                        example: `matrix m(1, 2);
m.arr = {{-0.0, 2.0}};
neg_zero(m);`,
                        notes: ["This is a formatting cleanup helper, not a tolerance cleanup helper."]
                    }),
                    createEntry({
                        id: "matrix-fpg",
                        title: "fpg",
                        signature: `void fpg(matrix& m)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Zeros matrix elements whose absolute value is at most EPS.",
                        example: `matrix m(1, 3);
m.arr = {{1e-8, -1e-7, 2.0}};
fpg(m);`,
                        notes: ["Used heavily after arithmetic and elimination to reduce printed residue."]
                    })
                ]
            },
            {
                title: "Elementary row and column operations",
                copy: "These helpers return transformed copies and leave the original matrix unchanged.",
                entries: [
                    createEntry({
                        id: "matrix-row-op",
                        title: "row_op",
                        signature: `matrix row_op(int i, double coeff, int j)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns a copy after applying row[i] += coeff * row[j].",
                        example: `matrix m(3, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};

matrix eliminated = m.row_op(1, -4.0, 0);`,
                        notes: ["Useful for explicit elimination steps."]
                    }),
                    createEntry({
                        id: "matrix-col-op",
                        title: "col_op",
                        signature: `matrix col_op(int i, double coeff, int j)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns a copy after applying col[i] += coeff * col[j].",
                        example: `matrix m(2, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

matrix changed = m.col_op(2, -1.0, 0);`,
                        notes: ["The original matrix remains unchanged."]
                    }),
                    createEntry({
                        id: "matrix-row-swap",
                        title: "row_swap",
                        signature: `matrix row_swap(int i, int j)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns a copy with rows i and j swapped.",
                        example: `matrix m(2, 2);
m.arr = {{1.0, 2.0}, {3.0, 4.0}};

matrix swapped = m.row_swap(0, 1);`,
                        notes: ["Row swaps are used internally by elimination routines when pivots are zero."]
                    }),
                    createEntry({
                        id: "matrix-col-swap",
                        title: "col_swap",
                        signature: `matrix col_swap(int i, int j)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns a copy with columns i and j swapped.",
                        example: `matrix m(2, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

matrix swapped = m.col_swap(0, 2);`,
                        notes: ["The returned copy is normalized with fpg()."]
                    }),
                    createEntry({
                        id: "matrix-row-multi",
                        title: "row_multi",
                        signature: `matrix row_multi(int i, double factor)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns a copy with row i scaled by factor.",
                        example: `matrix m(2, 2);
m.arr = {{1.0, 2.0}, {3.0, 4.0}};

matrix scaled = m.row_multi(1, 0.5);`,
                        notes: ["Common in Gaussian elimination and inverse calculation."]
                    }),
                    createEntry({
                        id: "matrix-col-multi",
                        title: "col_multi",
                        signature: `matrix col_multi(int i, double factor)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns a copy with column i scaled by factor.",
                        example: `matrix m(2, 2);
m.arr = {{1.0, 2.0}, {3.0, 4.0}};

matrix scaled = m.col_multi(0, -1.0);`,
                        notes: ["The input matrix itself is left untouched."]
                    })
                ]
            },
            {
                title: "Reduction routines",
                copy: "These methods mutate the matrix in place and report row-swap sign information.",
                entries: [
                    createEntry({
                        id: "matrix-echelon",
                        title: "echelon",
                        signature: `int echelon()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Transforms the matrix in place to row-echelon form and returns the sign from row swaps.",
                        example: `matrix m(3, 3);
m.arr = {{2.0, 1.0, -1.0}, {-3.0, -1.0, 2.0}, {-2.0, 1.0, 2.0}};

int sign = m.echelon();`,
                        notes: ["Pivot rows are not normalized to 1."]
                    }),
                    createEntry({
                        id: "matrix-gaussian",
                        title: "gaussian",
                        signature: `int gaussian()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Transforms the matrix in place with normalized pivots and clears entries below them.",
                        example: `matrix m(3, 4);
m.arr = {{2.0, 1.0, -1.0, 8.0}, {-3.0, -1.0, 2.0, -11.0}, {-2.0, 1.0, 2.0, -3.0}};

int sign = m.gaussian();`,
                        notes: ["Pivot rows are scaled so the pivot element becomes 1 when possible."]
                    }),
                    createEntry({
                        id: "matrix-gauss-jordan",
                        title: "gauss_jordan",
                        signature: `int gauss_jordan()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Applies Gaussian elimination and then clears entries above each pivot.",
                        example: `matrix m(3, 3);
m.arr = {{1.0, 2.0, -1.0}, {2.0, 4.0, 1.0}, {1.0, 0.0, 3.0}};

int sign = m.gauss_jordan();`,
                        notes: ["Useful when you want reduced row-echelon behavior."]
                    }),
                    createEntry({
                        id: "matrix-canonical",
                        title: "canonical",
                        signature: `int canonical()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Applies Gaussian elimination and then uses column operations to clear values to the right of pivots.",
                        example: `matrix m(2, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

int sign = m.canonical();`,
                        notes: ["This is a library-specific canonical-style reduction helper."]
                    }),
                    createEntry({
                        id: "matrix-rank",
                        title: "rank",
                        signature: `int rank()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Computes the rank by Gaussian elimination on a copy of the matrix.",
                        example: `matrix m(3, 3);
m.arr = {{1.0, 2.0, 3.0}, {2.0, 4.0, 6.0}, {1.0, 1.0, 1.0}};

int r = m.rank();`,
                        notes: ["The original matrix is not modified."]
                    }),
                    createEntry({
                        id: "matrix-free-variable",
                        title: "free_variable",
                        signature: `int free_variable()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns row - rank() for a square coefficient matrix.",
                        example: `matrix coeff(3, 3);
coeff.arr = {{1.0, 2.0, 3.0}, {2.0, 4.0, 6.0}, {0.0, 0.0, 0.0}};

int free_vars = coeff.free_variable();`,
                        notes: ["Throws when the matrix is not square."]
                    })
                ]
            },
            {
                title: "Properties, derived matrices, and solving",
                copy: "These APIs compute matrix invariants, related matrices, and linear-system solutions.",
                entries: [
                    createEntry({
                        id: "matrix-trace",
                        title: "trace",
                        signature: `double trace()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the sum of the main diagonal of a square matrix.",
                        example: `matrix m(3, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};

double tr = m.trace();`,
                        notes: ["Throws when the matrix is not square."]
                    }),
                    createEntry({
                        id: "matrix-diag-product",
                        title: "diag_product",
                        signature: `double diag_product()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the product of the main diagonal entries over min(row, col).",
                        example: `matrix m(3, 3);
m.arr = {{2.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {0.0, 0.0, 4.0}};

double prod = m.diag_product();`,
                        notes: ["Often used after echelon reduction while computing determinants."]
                    }),
                    createEntry({
                        id: "matrix-det",
                        title: "det",
                        signature: `double det()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the determinant using echelon reduction and diagonal multiplication.",
                        example: `matrix m(3, 3);
m.arr = {{1.0, 2.0, 3.0}, {0.0, 1.0, 4.0}, {5.0, 6.0, 0.0}};

double determinant = m.det();`,
                        notes: ["Defined only for square matrices."]
                    }),
                    createEntry({
                        id: "matrix-minor",
                        title: "minor",
                        signature: `double minor(int x, int y)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Computes the minor of entry (x, y) by removing its row and column.",
                        example: `matrix m(3, 3);
m.arr = {{3.0, 0.0, 2.0}, {2.0, 0.0, -2.0}, {0.0, 1.0, 1.0}};

double mn = m.minor(0, 0);`,
                        notes: ["Square matrices only."]
                    }),
                    createEntry({
                        id: "matrix-cofactor",
                        title: "cofactor",
                        signature: `double cofactor(int x, int y)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the signed minor corresponding to entry (x, y).",
                        example: `matrix m(3, 3);
m.arr = {{3.0, 0.0, 2.0}, {2.0, 0.0, -2.0}, {0.0, 1.0, 1.0}};

double cof = m.cofactor(0, 1);`,
                        notes: ["The sign follows the usual checkerboard pattern."]
                    }),
                    createEntry({
                        id: "matrix-adjoint",
                        title: "adjoint",
                        signature: `matrix adjoint()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Builds the classical adjoint by computing cofactors and transposing them.",
                        example: `matrix m(3, 3);
m.arr = {{4.0, 7.0, 2.0}, {3.0, 6.0, 1.0}, {2.0, 5.0, 1.0}};

matrix adj = m.adjoint();`,
                        notes: ["Defined only for square matrices."]
                    }),
                    createEntry({
                        id: "matrix-transpose",
                        title: "transpose",
                        signature: `matrix transpose()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the transposed matrix.",
                        example: `matrix m(2, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

matrix t = m.transpose();`,
                        notes: ["The original matrix is not modified."]
                    }),
                    createEntry({
                        id: "matrix-inverse",
                        title: "inverse",
                        signature: `matrix inverse()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the inverse of a nonsingular square matrix using Gauss-Jordan elimination.",
                        example: `matrix m(3, 3);
m.arr = {{1.0, 2.0, 3.0}, {0.0, 1.0, 4.0}, {5.0, 6.0, 0.0}};

matrix inv = m.inverse();`,
                        notes: ["Throws if the matrix is singular or not square."]
                    }),
                    createEntry({
                        id: "matrix-solve",
                        title: "solve",
                        signature: `matrix solve()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Solves an N x (N + 1) augmented system and returns the solution column vector.",
                        example: `matrix aug(3, 4);
aug.arr = {
    {2.0, 1.0, -1.0, 8.0},
    {-3.0, -1.0, 2.0, -11.0},
    {-2.0, 1.0, 2.0, -3.0}
};

matrix solution = aug.solve();`,
                        notes: ["The implementation assumes a solvable augmented system with a usable pivot structure."]
                    }),
                    createEntry({
                        id: "matrix-norm",
                        title: "norm",
                        signature: `double norm()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the Frobenius norm of the matrix.",
                        example: `matrix m(2, 2);
m.arr = {{1.0, 2.0}, {3.0, 4.0}};

double n = m.norm();`,
                        notes: ["Internally implemented as sqrt(inner_product(*this, *this))."]
                    })
                ]
            },
            {
                title: "Views and replacements",
                copy: "These helpers extract or overwrite structural slices.",
                entries: [
                    createEntry({
                        id: "matrix-get-row",
                        title: "get_row",
                        signature: `matrix get_row(int r)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns a same-sized matrix whose row r matches the input and whose other rows are zero.",
                        example: `matrix m(3, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};

matrix only_second_row = m.get_row(1);`,
                        notes: ["This is not the same as get_row_vec()."]
                    }),
                    createEntry({
                        id: "matrix-get-col",
                        title: "get_col",
                        signature: `matrix get_col(int c)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns a same-sized matrix whose column c matches the input and whose other columns are zero.",
                        example: `matrix m(3, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};

matrix only_first_col = m.get_col(0);`,
                        notes: ["Useful inside orthogonalization routines."]
                    }),
                    createEntry({
                        id: "matrix-get-row-vec",
                        title: "get_row_vec",
                        signature: `matrix get_row_vec(int r)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns row r as a 1 x col matrix.",
                        example: `matrix m(2, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

matrix row_vec = m.get_row_vec(0);`,
                        notes: ["This is a compact row vector representation."]
                    }),
                    createEntry({
                        id: "matrix-get-col-vec",
                        title: "get_col_vec",
                        signature: `matrix get_col_vec(int c)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns column c as a row x 1 matrix.",
                        example: `matrix m(2, 3);
m.arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

matrix col_vec = m.get_col_vec(2);`,
                        notes: ["This is the column-vector form used by the Gram-Schmidt helpers."]
                    }),
                    createEntry({
                        id: "matrix-replace-row",
                        title: "replace_row",
                        signature: `void replace_row(int r, matrix rw)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Copies the first row of rw into row r of the current matrix.",
                        example: `matrix m(2, 3);
matrix replacement(1, 3);
replacement.arr = {{9.0, 9.0, 9.0}};

m.replace_row(0, replacement);`,
                        notes: ["The caller is responsible for providing a compatible source row."]
                    }),
                    createEntry({
                        id: "matrix-replace-col",
                        title: "replace_col",
                        signature: `void replace_col(int c, matrix cn)`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Copies the first column of cn into column c of the current matrix.",
                        example: `matrix m(3, 2);
matrix replacement(3, 1);
replacement.arr = {{7.0}, {8.0}, {9.0}};

m.replace_col(1, replacement);`,
                        notes: ["The source matrix is expected to expose a usable first column."]
                    })
                ]
            },
            {
                title: "Decompositions and spectral helpers",
                copy: "These methods support orthogonalization, factorization, characteristic polynomials, and eigenvalue estimation.",
                entries: [
                    createEntry({
                        id: "matrix-orthogonalize",
                        title: "orthogonalize",
                        signature: `matrix orthogonalize()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Applies classical Gram-Schmidt to the columns and returns orthogonal output columns.",
                        example: `matrix basis(3, 2);
basis.arr = {{1.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}};

matrix ortho = basis.orthogonalize();`,
                        notes: ["The returned matrix has the same shape as the input."]
                    }),
                    createEntry({
                        id: "matrix-orthonormalize",
                        title: "orthonormalize",
                        signature: `matrix orthonormalize()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Applies classical Gram-Schmidt and normalizes each resulting column.",
                        example: `matrix basis(3, 2);
basis.arr = {{1.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}};

matrix q = basis.orthonormalize();`,
                        notes: ["Equivalent to the Q factor used by qr_decomp_q()."]
                    }),
                    createEntry({
                        id: "matrix-qr-q",
                        title: "qr_decomp_q",
                        signature: `matrix qr_decomp_q()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the orthonormal Q factor from a QR decomposition.",
                        example: `matrix a(3, 3);
a.arr = {{12.0, -51.0, 4.0}, {6.0, 167.0, -68.0}, {-4.0, 24.0, -41.0}};

matrix q = a.qr_decomp_q();`,
                        notes: ["Built from orthonormalize()."]
                    }),
                    createEntry({
                        id: "matrix-qr-r",
                        title: "qr_decomp_r",
                        signature: `matrix qr_decomp_r()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the upper-triangular R factor from a QR decomposition.",
                        example: `matrix a(3, 3);
a.arr = {{12.0, -51.0, 4.0}, {6.0, 167.0, -68.0}, {-4.0, 24.0, -41.0}};

matrix r = a.qr_decomp_r();`,
                        notes: ["Computed from the Q factor and the original columns."]
                    }),
                    createEntry({
                        id: "matrix-lu-decomposition",
                        title: "lu_decomposition",
                        signature: `std::tuple<matrix, matrix, matrix> lu_decomposition()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Performs LU decomposition with partial pivoting: P * A = L * U.",
                        example: `matrix a(3, 3);
a.arr = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};
auto [p, l, u] = a.lu_decomposition();`,
                        notes: ["Returns a tuple of {P, L, U}.", "Throws std::runtime_error if the matrix is singular or rectangular."]
                    }),
                    createEntry({
                        id: "matrix-char-poly",
                        title: "char_poly",
                        signature: `matrix char_poly()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns characteristic polynomial coefficients as a 1 x (n + 1) row vector.",
                        example: `matrix a(3, 3);
a.arr = {{1.0, 2.0, 3.0}, {0.0, 1.0, 4.0}, {5.0, 6.0, 0.0}};

matrix coeffs = a.char_poly();`,
                        notes: ["The implementation uses a Leverrier-Faddeev style recurrence."]
                    }),
                    createEntry({
                        id: "matrix-eigenvalues",
                        title: "eigenvalues",
                        signature: `matrix eigenvalues()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Estimates real eigenvalues using repeated unshifted QR iteration and returns them as a column vector.",
                        example: `matrix a(3, 3);
a.arr = {{3.0, 4.0, 5.0}, {6.0, 7.0, 8.0}, {8.0, 2.0, 3.0}};

matrix lambda = a.eigenvalues();`,
                        notes: ["Iteration stops when the working matrix becomes exactly upper triangular after cleanup."]
                    }),
                    createEntry({
                        id: "matrix-upper-tri-check",
                        title: "check_upper_tri",
                        signature: `bool check_upper_tri()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns true when all strictly lower-triangular entries are exactly zero.",
                        example: `matrix u(3, 3);
u.arr = {{1.0, 2.0, 3.0}, {0.0, 4.0, 5.0}, {0.0, 0.0, 6.0}};

bool is_upper = u.check_upper_tri();`,
                        notes: ["This is an exact check rather than a tolerance-based one."]
                    }),
                    createEntry({
                        id: "matrix-lower-tri-check",
                        title: "check_lower_tri",
                        signature: `bool check_lower_tri()`,
                        kind: "method",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns true when all strictly upper-triangular entries are exactly zero.",
                        example: `matrix l(3, 3);
l.arr = {{1.0, 0.0, 0.0}, {2.0, 3.0, 0.0}, {4.0, 5.0, 6.0}};

bool is_lower = l.check_lower_tri();`,
                        notes: ["This is also an exact check."]
                    })
                ]
            },
            {
                title: "Free matrix algorithms",
                copy: "These helpers operate on matrix objects outside the class itself.",
                entries: [
                    createEntry({
                        id: "matrix-matpow",
                        title: "matpow",
                        signature: `matrix matpow(matrix mat, long long expo)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Raises a square matrix to a non-negative integer power using repeated squaring.",
                        example: `matrix a(2, 2);
a.arr = {{1.0, 1.0}, {1.0, 0.0}};

matrix a5 = matpow(a, 5);`,
                        notes: [
                            "Throws std::runtime_error for non-square matrices and negative exponents."
                        ]
                    }),
                    createEntry({
                        id: "matrix-check-ortho",
                        title: "check_ortho",
                        signature: `bool check_ortho(matrix mat)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Checks whether transpose() and inverse() are exactly equal.",
                        example: `matrix q(2, 2);
q.arr = {{0.0, -1.0}, {1.0, 0.0}};

bool orthogonal = check_ortho(q);`,
                        notes: ["Sensitive to floating-point roundoff because equality is exact."]
                    }),
                    createEntry({
                        id: "matrix-check-unitary",
                        title: "check_unitary",
                        signature: `bool check_unitary(matrix mat)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Legacy helper that checks whether transpose() and adjoint() are exactly equal.",
                        example: `matrix m = identity(3);
bool unitary_like = check_unitary(m);`,
                        notes: ["This helper uses the classical adjoint, not a complex conjugate transpose."]
                    }),
                    createEntry({
                        id: "matrix-inner-product",
                        title: "inner_product",
                        signature: `double inner_product(matrix a, matrix b)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Computes the Frobenius inner product of two same-shaped matrices.",
                        example: `matrix a(2, 1);
matrix b(2, 1);
a.arr = {{1.0}, {2.0}};
b.arr = {{3.0}, {4.0}};

double dot = inner_product(a, b);`,
                        notes: ["Throws when the shapes differ."]
                    }),
                    createEntry({
                        id: "matrix-angle",
                        title: "angle",
                        signature: `double angle(matrix a, matrix b)`,
                        kind: "free function",
                        declaredIn: "include/elda/matrix.hpp",
                        implementedIn: "src/matrix.cpp",
                        description: "Returns the angle in radians between two matrices under the Frobenius inner product.",
                        example: `matrix a(2, 1);
matrix b(2, 1);
a.arr = {{1.0}, {0.0}};
b.arr = {{0.0}, {1.0}};

double theta = angle(a, b);`,
                        notes: ["Equivalent to the familiar vector-angle formula in Frobenius space."]
                    })
                ]
            }
        ]
    },
    {
        id: "transforms",
        mountId: "api-transforms",
        title: "Homogeneous transforms",
        headerPath: "include/elda/transforms.hpp",
        synopsis: `#include <elda/transforms.hpp>`,
        intro:
            "This header declares the 2D and 3D homogeneous transformation builders returned as matrix objects.",
        chips: ["2D transforms", "3D transforms", "radians"],
        groups: [
            {
                title: "2D homogeneous transforms",
                copy: "These functions return 3 x 3 matrices suitable for homogeneous 2D coordinates.",
                entries: [
                    createEntry({
                        id: "transforms-shift-2d",
                        title: "shift (2D)",
                        signature: `matrix shift(double dx, double dy)`,
                        kind: "free function",
                        declaredIn: "include/elda/transforms.hpp",
                        implementedIn: "src/transforms.cpp",
                        description: "Builds a 3 x 3 translation matrix that shifts x by dx and y by dy.",
                        example: `matrix t = shift(2.0, -1.0);
matrix point = vec3(1.0, 2.0, 1.0);
matrix moved = t * point;`,
                        notes: ["The translation is stored in the last column."]
                    }),
                    createEntry({
                        id: "transforms-scale-2d",
                        title: "scale (2D)",
                        signature: `matrix scale(double kx, double ky)`,
                        kind: "free function",
                        declaredIn: "include/elda/transforms.hpp",
                        implementedIn: "src/transforms.cpp",
                        description: "Builds a 3 x 3 scaling matrix for homogeneous 2D coordinates.",
                        example: `matrix s = scale(2.0, 0.5);
matrix point = vec3(3.0, 4.0, 1.0);
matrix scaled = s * point;`,
                        notes: ["Scale factors occupy the upper-left diagonal entries."]
                    }),
                    createEntry({
                        id: "transforms-rotate-2d",
                        title: "rotate",
                        signature: `matrix rotate(double angle)`,
                        kind: "free function",
                        declaredIn: "include/elda/transforms.hpp",
                        implementedIn: "src/transforms.cpp",
                        description: "Builds a 3 x 3 rotation matrix for homogeneous 2D coordinates.",
                        example: `matrix r = rotate(PI / 2.0);
matrix point = vec3(1.0, 0.0, 1.0);
matrix turned = r * point;`,
                        notes: ["Angles are interpreted in radians."]
                    })
                ]
            },
            {
                title: "3D homogeneous transforms",
                copy: "These functions return 4 x 4 matrices suitable for homogeneous 3D coordinates.",
                entries: [
                    createEntry({
                        id: "transforms-shift-3d",
                        title: "shift (3D)",
                        signature: `matrix shift(double dx, double dy, double dz)`,
                        kind: "free function",
                        declaredIn: "include/elda/transforms.hpp",
                        implementedIn: "src/transforms.cpp",
                        description: "Builds a 4 x 4 translation matrix for homogeneous 3D coordinates.",
                        example: `matrix t = shift(1.0, 2.0, 3.0);`,
                        notes: ["Translation terms are placed in the last column."]
                    }),
                    createEntry({
                        id: "transforms-scale-3d",
                        title: "scale (3D)",
                        signature: `matrix scale(double kx, double ky, double kz)`,
                        kind: "free function",
                        declaredIn: "include/elda/transforms.hpp",
                        implementedIn: "src/transforms.cpp",
                        description: "Builds a 4 x 4 scaling matrix for homogeneous 3D coordinates.",
                        example: `matrix s = scale(2.0, 2.0, 0.5);`,
                        notes: ["The scale factors occupy the x, y, and z diagonal positions."]
                    }),
                    createEntry({
                        id: "transforms-rot-x",
                        title: "rot_x",
                        signature: `matrix rot_x(double angle)`,
                        kind: "free function",
                        declaredIn: "include/elda/transforms.hpp",
                        implementedIn: "src/transforms.cpp",
                        description: "Builds a 4 x 4 rotation matrix around the x-axis.",
                        example: `matrix rx = rot_x(PI / 6.0);`,
                        notes: ["Angle values are measured in radians."]
                    }),
                    createEntry({
                        id: "transforms-rot-y",
                        title: "rot_y",
                        signature: `matrix rot_y(double angle)`,
                        kind: "free function",
                        declaredIn: "include/elda/transforms.hpp",
                        implementedIn: "src/transforms.cpp",
                        description: "Builds a 4 x 4 rotation matrix around the y-axis.",
                        example: `matrix ry = rot_y(PI / 4.0);`,
                        notes: ["The upper-left and lower-left blocks carry cos/sin terms."]
                    }),
                    createEntry({
                        id: "transforms-rot-z",
                        title: "rot_z",
                        signature: `matrix rot_z(double angle)`,
                        kind: "free function",
                        declaredIn: "include/elda/transforms.hpp",
                        implementedIn: "src/transforms.cpp",
                        description: "Builds a 4 x 4 rotation matrix around the z-axis.",
                        example: `matrix rz = rot_z(PI / 3.0);`,
                        notes: ["Equivalent in structure to the 2D rotation helper, embedded in 3D homogeneous form."]
                    })
                ]
            }
        ]
    },
    {
        id: "vector-utils",
        mountId: "api-vector-utils",
        title: "Vector utility helpers",
        headerPath: "include/elda/vector_utils.hpp",
        synopsis: `#include <elda/vector_utils.hpp>`,
        intro:
            "This header provides small column-vector factory helpers and linear-combination checks built on rank comparisons.",
        chips: ["column vectors", "span checks"],
        groups: [
            {
                title: "Zero-initialized vector factories",
                copy: "These helpers return small zero column vectors without manual shape construction.",
                entries: [
                    createEntry({
                        id: "vector-vec1-zero",
                        title: "vec1",
                        signature: `matrix vec1()`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a zero-initialized 1 x 1 column vector.",
                        example: `matrix v = vec1();`,
                        notes: ["Equivalent to matrix(1, 1)."]
                    }),
                    createEntry({
                        id: "vector-vec2-zero",
                        title: "vec2",
                        signature: `matrix vec2()`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a zero-initialized 2 x 1 column vector.",
                        example: `matrix v = vec2();`,
                        notes: ["Equivalent to matrix(2, 1)."]
                    }),
                    createEntry({
                        id: "vector-vec3-zero",
                        title: "vec3",
                        signature: `matrix vec3()`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a zero-initialized 3 x 1 column vector.",
                        example: `matrix v = vec3();`,
                        notes: ["Useful for homogeneous 2D points before assigning coordinates."]
                    }),
                    createEntry({
                        id: "vector-vec4-zero",
                        title: "vec4",
                        signature: `matrix vec4()`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a zero-initialized 4 x 1 column vector.",
                        example: `matrix v = vec4();`,
                        notes: ["Useful for homogeneous 3D coordinates."]
                    }),
                    createEntry({
                        id: "vector-vec5-zero",
                        title: "vec5",
                        signature: `matrix vec5()`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a zero-initialized 5 x 1 column vector.",
                        example: `matrix v = vec5();`,
                        notes: ["Equivalent to matrix(5, 1)."]
                    })
                ]
            },
            {
                title: "Coordinate vector factories",
                copy: "These overloads return populated column vectors from explicit coordinates.",
                entries: [
                    createEntry({
                        id: "vector-vec1-values",
                        title: "vec1(x)",
                        signature: `matrix vec1(double x)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a 1 x 1 column vector containing x.",
                        example: `matrix v = vec1(5.0);`,
                        notes: ["The single coordinate is stored in arr[0][0]."]
                    }),
                    createEntry({
                        id: "vector-vec2-values",
                        title: "vec2(x, y)",
                        signature: `matrix vec2(double x, double y)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a 2 x 1 column vector containing x and y.",
                        example: `matrix v = vec2(3.0, 4.0);`,
                        notes: ["Common for planar vector examples."]
                    }),
                    createEntry({
                        id: "vector-vec3-values",
                        title: "vec3(x, y, z)",
                        signature: `matrix vec3(double x, double y, double z)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a 3 x 1 column vector containing x, y, and z.",
                        example: `matrix point = vec3(1.0, 2.0, 1.0);`,
                        notes: ["Often used as a homogeneous 2D point with final entry 1."]
                    }),
                    createEntry({
                        id: "vector-vec4-values",
                        title: "vec4(x, y, z, w)",
                        signature: `matrix vec4(double x, double y, double z, double w)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a 4 x 1 column vector containing x, y, z, and w.",
                        example: `matrix point = vec4(1.0, 2.0, 3.0, 1.0);`,
                        notes: ["Suitable for homogeneous 3D points."]
                    }),
                    createEntry({
                        id: "vector-vec5-values",
                        title: "vec5(x, y, z, w, u)",
                        signature: `matrix vec5(double x, double y, double z, double w, double u)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns a 5 x 1 column vector containing five coordinates.",
                        example: `matrix v = vec5(1.0, 2.0, 3.0, 4.0, 5.0);`,
                        notes: ["Useful for small higher-dimensional examples."]
                    })
                ]
            },
            {
                title: "Linear-combination checks",
                copy: "These overloads test whether the final vector belongs to the span of the earlier vectors.",
                entries: [
                    createEntry({
                        id: "vector-check-lin-comb-2",
                        title: "check_lin_comb(m1, m2)",
                        signature: `bool check_lin_comb(matrix m1, matrix m2)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns true when m2 lies in the span of m1.",
                        example: `matrix v1 = vec2(1.0, 2.0);
matrix v2 = vec2(2.0, 4.0);

bool in_span = check_lin_comb(v1, v2);`,
                        notes: ["Implemented by comparing the ranks of [m1 m2] and [m1]."]
                    }),
                    createEntry({
                        id: "vector-check-lin-comb-3",
                        title: "check_lin_comb(m1, m2, m3)",
                        signature: `bool check_lin_comb(matrix m1, matrix m2, matrix m3)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns true when m3 lies in the span of m1 and m2.",
                        example: `matrix v1 = vec3(1.0, 0.0, 0.0);
matrix v2 = vec3(0.0, 1.0, 0.0);
matrix v3 = vec3(2.0, 3.0, 0.0);

bool in_span = check_lin_comb(v1, v2, v3);`,
                        notes: ["The candidate vector is always the last argument."]
                    }),
                    createEntry({
                        id: "vector-check-lin-comb-4",
                        title: "check_lin_comb(m1, m2, m3, m4)",
                        signature: `bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns true when m4 lies in the span of m1, m2, and m3.",
                        example: `matrix v1 = vec3(1.0, 0.0, 0.0);
matrix v2 = vec3(0.0, 1.0, 0.0);
matrix v3 = vec3(0.0, 0.0, 1.0);
matrix v4 = vec3(2.0, 3.0, 4.0);

bool in_span = check_lin_comb(v1, v2, v3, v4);`,
                        notes: ["The implementation constructs an augmented matrix of the input columns."]
                    }),
                    createEntry({
                        id: "vector-check-lin-comb-5",
                        title: "check_lin_comb(m1, m2, m3, m4, m5)",
                        signature: `bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4, matrix m5)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns true when m5 lies in the span of m1, m2, m3, and m4.",
                        example: `matrix v1 = vec4(1.0, 0.0, 0.0, 0.0);
matrix v2 = vec4(0.0, 1.0, 0.0, 0.0);
matrix v3 = vec4(0.0, 0.0, 1.0, 0.0);
matrix v4 = vec4(0.0, 0.0, 0.0, 1.0);
matrix v5 = vec4(1.0, 2.0, 3.0, 4.0);

bool in_span = check_lin_comb(v1, v2, v3, v4, v5);`,
                        notes: ["With a full basis, the final vector belongs to the span by construction."]
                    }),
                    createEntry({
                        id: "vector-check-lin-comb-list",
                        title: "check_lin_comb(target, vectors)",
                        signature: `bool check_lin_comb(matrix target, const std::vector<matrix>& vectors)`,
                        kind: "free function",
                        declaredIn: "include/elda/vector_utils.hpp",
                        implementedIn: "src/vector_utils.cpp",
                        description: "Returns true when target lies in the span of the given vectors.",
                        example: `matrix v1 = vec3(1.0, 0.0, 0.0);
matrix v2 = vec3(0.0, 1.0, 0.0);
matrix target = vec3(5.0, 5.0, 0.0);
std::vector<matrix> vec_set = {v1, v2};

bool in_span = check_lin_comb(target, vec_set);`,
                        notes: ["This overload handles an arbitrary number of input vectors."]
                    })
                ]
            }
        ]
    }
];

const exampleData = {
    arithmetic: {
        title: "3x3 arithmetic",
        description:
            "Matrix addition, multiplication, and transpose are direct copy-returning operations on linalg::matrix. The values below are produced from the current ELDA build.",
        bullets: [
            "Addition and subtraction require identical shapes.",
            "Matrix multiplication requires lhs.col == rhs.row.",
            "Arithmetic results are normalized with fpg() before return."
        ],
        code: `matrix a(3, 3);
a.arr = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};

matrix b(3, 3);
b.arr = {{-2, 1, 0}, {3, 0, 1}, {4, -1, 2}};

matrix sum = a + b;
matrix product = a * b;
matrix transposed = a.transpose();`,
        visuals: [
            {
                type: "equation",
                heading: "Addition",
                lhs: [[1, 2, 3], [0, 1, 4], [5, 6, 0]],
                op: "+",
                rhs: [[-2, 1, 0], [3, 0, 1], [4, -1, 2]],
                result: [[-1, 3, 3], [3, 1, 5], [9, 5, 2]],
                caption: "A + B"
            },
            {
                type: "equation",
                heading: "Multiplication",
                lhs: [[1, 2, 3], [0, 1, 4], [5, 6, 0]],
                op: "x",
                rhs: [[-2, 1, 0], [3, 0, 1], [4, -1, 2]],
                result: [[16, -2, 8], [19, -4, 9], [8, 5, 6]],
                caption: "A x B"
            },
            {
                type: "single",
                heading: "Transpose",
                matrix: [[1, 0, 5], [2, 1, 6], [3, 4, 0]],
                caption: "transpose(A)"
            }
        ]
    },
    inverse: {
        title: "Determinant and inverse",
        description:
            "For square matrices, ELDA provides determinant, adjoint, inverse, characteristic polynomial, and eigenvalue helpers.",
        bullets: [
            "det() works through echelon reduction and diagonal multiplication.",
            "inverse() uses Gauss-Jordan elimination on [A | I].",
            "The input matrix must be square and nonsingular."
        ],
        code: `matrix a(3, 3);
a.arr = {{1, 2, 3}, {0, 1, 4}, {5, 6, 0}};

double determinant = a.det();
matrix inverse = a.inverse();
matrix coefficients = a.char_poly();
matrix lambda = a.eigenvalues();`,
        visuals: [
            {
                type: "single",
                heading: "Input matrix",
                matrix: [[1, 2, 3], [0, 1, 4], [5, 6, 0]],
                caption: "A, with det(A) = 1"
            },
            {
                type: "single",
                heading: "Inverse",
                matrix: [[-24, 18, 5], [20, -15, -4], [-5, 4, 1]],
                caption: "A^-1 returned by matrix::inverse()"
            },
            {
                type: "single",
                heading: "Characteristic polynomial coefficients",
                matrix: [[1, -2, -38, -1]],
                caption: "char_poly(A)"
            },
            {
                type: "single",
                heading: "Eigenvalue estimates",
                matrix: [[7.256021], [-5.229669], [-0.026323]],
                caption: "eigenvalues(A)"
            }
        ]
    },
    solve: {
        title: "Solving a 3x3 linear system",
        description:
            "solve() expects an augmented matrix with shape N x (N + 1). The coefficient matrix is reduced with gaussian(), then the solution vector is recovered by back-substitution.",
        bullets: [
            "The returned object is an N x 1 column vector.",
            "The implementation assumes a solvable augmented system.",
            "This example produces x = [2, 3, -1]^T."
        ],
        code: `matrix aug(3, 4);
aug.arr = {
    { 2,  1, -1,  8},
    {-3, -1,  2, -11},
    {-2,  1,  2, -3}
};

matrix solution = aug.solve();`,
        visuals: [
            {
                type: "equation",
                heading: "Coefficient system",
                lhs: [[2, 1, -1], [-3, -1, 2], [-2, 1, 2]],
                op: "x",
                rhs: [[2], [3], [-1]],
                result: [[8], [-11], [-3]],
                caption: "A x = b"
            },
            {
                type: "single",
                heading: "Recovered solution vector",
                matrix: [[2], [3], [-1]],
                caption: "solve(aug)"
            }
        ]
    },
    qr: {
        title: "QR decomposition on a 3x3 matrix",
        description:
            "ELDA builds QR through classical Gram-Schmidt orthonormalization of columns.",
        bullets: [
            "qr_decomp_q() returns the orthonormal factor.",
            "qr_decomp_r() returns the upper-triangular factor.",
            "The matrix below is a standard QR benchmark."
        ],
        code: `matrix a(3, 3);
a.arr = {{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}};

matrix q = a.qr_decomp_q();
matrix r = a.qr_decomp_r();`,
        visuals: [
            {
                type: "single",
                heading: "Q factor",
                matrix: [
                    [0.857143, -0.394286, -0.331429],
                    [0.428571, 0.902857, 0.034286],
                    [-0.285714, 0.171429, -0.942857]
                ],
                caption: "qr_decomp_q()"
            },
            {
                type: "single",
                heading: "R factor",
                matrix: [[14, 21, -14], [0, 175, -70], [0, 0, 35]],
                caption: "qr_decomp_r()"
            }
        ]
    },
    transform: {
        title: "3x3 homogeneous 2D transforms",
        description:
            "In 2D, translation and rotation are represented as 3 x 3 homogeneous matrices, making composition straightforward.",
        bullets: [
            "shift(dx, dy) writes the translation into the last column.",
            "rotate(angle) stores cos/sin terms in the upper-left 2 x 2 block.",
            "These matrices are ready for homogeneous point vectors [x, y, 1]^T."
        ],
        code: `matrix translation = shift(2.0, -1.0);
matrix quarter_turn = rotate(PI / 2.0);

matrix point = vec3(1.0, 2.0, 1.0);
matrix moved = translation * point;`,
        visuals: [
            {
                type: "single",
                heading: "Translation matrix",
                matrix: [[1, 0, 2], [0, 1, -1], [0, 0, 1]],
                caption: "shift(2, -1)"
            },
            {
                type: "single",
                heading: "Rotation matrix",
                matrix: [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                caption: "Mathematical form of rotate(PI / 2); raw printed output may show signed zeros"
            },
            {
                type: "equation",
                heading: "Point translation",
                lhs: [[1, 0, 2], [0, 1, -1], [0, 0, 1]],
                op: "x",
                rhs: [[1], [2], [1]],
                result: [[3], [1], [1]],
                caption: "Applying shift(2, -1) to [1, 2, 1]^T"
            }
        ]
    }
};

function escapeHtml(value) {
    return value
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;");
}

function renderHeaderSection(header) {
    const groupsHtml = header.groups.map(renderGroup).join("");
    return `
        <article class="header-card">
            <h4>${header.title}</h4>
            <p>${header.intro}</p>
            <pre><code class="language-cpp">${escapeHtml(header.synopsis)}</code></pre>
            <div class="header-meta">
                <span class="chip">${header.groups.reduce((count, group) => count + group.entries.length, 0)} entries</span>
                ${header.chips.map((chip) => `<span class="chip">${chip}</span>`).join("")}
            </div>
        </article>
        <div class="api-groups">${groupsHtml}</div>
    `;
}

function renderGroup(group) {
    const entriesHtml = group.entries.map(renderEntry).join("");
    return `
        <section class="api-group" data-group-title="${escapeHtml(group.title)}">
            <h5>${group.title}</h5>
            <p class="api-group-copy">${group.copy}</p>
            <div class="api-card-stack">
                ${entriesHtml}
            </div>
        </section>
    `;
}

function renderEntry(entry) {
    const notesHtml = entry.notes.length
        ? `<ul>${entry.notes.map((note) => `<li>${note}</li>`).join("")}</ul>`
        : `<p>No additional notes.</p>`;

    const searchText = [
        entry.title,
        entry.signature,
        entry.description,
        entry.declaredIn,
        entry.implementedIn,
        entry.example,
        entry.notes.join(" ")
    ]
        .join(" ")
        .toLowerCase();

    return `
        <details class="api-card" data-search="${escapeHtml(searchText)}">
            <summary>
                <div>
                    <h6 class="api-card-title">${entry.title}</h6>
                    <p class="api-card-signature"><code>${escapeHtml(entry.signature)}</code></p>
                </div>
                <div class="api-card-meta">
                    <span class="chip">${entry.kind}</span>
                </div>
            </summary>
            <div class="api-card-body">
                <div class="api-body-grid">
                    <div class="api-panel">
                        <h6>Declaration</h6>
                        <pre><code class="language-cpp">${escapeHtml(entry.signature)}</code></pre>
                        <h6>Description</h6>
                        <p>${entry.description}</p>
                        <div class="api-meta-list">
                            <div><strong>Declared in:</strong> <code>${entry.declaredIn}</code></div>
                            <div><strong>Implemented in:</strong> <code>${entry.implementedIn}</code></div>
                        </div>
                    </div>
                    <div class="api-panel">
                        <h6>Example</h6>
                        <pre><code class="language-cpp">${escapeHtml(entry.example)}</code></pre>
                        <h6>Notes</h6>
                        ${notesHtml}
                    </div>
                </div>
            </div>
        </details>
    `;
}

function renderReference() {
    referenceData.forEach((header) => {
        const mount = document.getElementById(header.mountId);
        if (!mount) {
            return;
        }
        mount.innerHTML = renderHeaderSection(header);
    });
}

function formatValue(value) {
    if (Number.isInteger(value)) {
        return String(value);
    }
    const normalized = Number(value.toFixed(6));
    return String(normalized);
}

function renderMatrix(matrix) {
    const rows = matrix.length;
    const cols = matrix[0].length;
    const cells = matrix
        .flat()
        .map((value) => `<span>${formatValue(value)}</span>`)
        .join("");

    return `
        <figure class="matrix-figure" aria-label="Matrix with ${rows} rows and ${cols} columns">
            <div class="matrix-grid" style="grid-template-columns: repeat(${cols}, minmax(2.4rem, auto));">
                ${cells}
            </div>
        </figure>
    `;
}

function renderVisual(visual) {
    if (visual.type === "equation") {
        return `
            <section class="visual-block">
                <h5>${visual.heading}</h5>
                <div class="equation-row">
                    ${renderMatrix(visual.lhs)}
                    <span class="equation-symbol">${visual.op}</span>
                    ${renderMatrix(visual.rhs)}
                    <span class="equation-symbol">=</span>
                    ${renderMatrix(visual.result)}
                </div>
                <p class="visual-caption">${visual.caption}</p>
            </section>
        `;
    }

    return `
        <section class="visual-block">
            <h5>${visual.heading}</h5>
            ${renderMatrix(visual.matrix)}
            <p class="visual-caption">${visual.caption}</p>
        </section>
    `;
}

function activateExample(exampleId) {
    const example = exampleData[exampleId];
    if (!example) {
        return;
    }

    document.getElementById("example-title").textContent = example.title;
    document.getElementById("example-description").textContent = example.description;
    document.getElementById("example-code").textContent = example.code;
    document.getElementById("example-bullets").innerHTML = example.bullets
        .map((bullet) => `<li>${bullet}</li>`)
        .join("");
    document.getElementById("example-visuals").innerHTML = example.visuals.map(renderVisual).join("");

    document.querySelectorAll(".example-tab").forEach((button) => {
        const active = button.dataset.example === exampleId;
        button.classList.toggle("active", active);
        button.setAttribute("aria-selected", active ? "true" : "false");
    });
}

function setupExampleTabs() {
    document.querySelectorAll(".example-tab").forEach((button) => {
        button.addEventListener("click", () => activateExample(button.dataset.example));
    });

    activateExample("arithmetic");
}

function setupApiSearch() {
    const search = document.getElementById("api-search");
    if (!search) {
        return;
    }

    search.addEventListener("input", () => {
        const query = search.value.trim().toLowerCase();

        document.querySelectorAll(".api-card").forEach((card) => {
            const match = card.dataset.search.includes(query);
            card.hidden = !match;
            if (query && match) {
                card.open = true;
            }
        });

        document.querySelectorAll(".api-group").forEach((group) => {
            const hasVisibleCards = Array.from(group.querySelectorAll(".api-card")).some((card) => !card.hidden);
            group.hidden = !hasVisibleCards;
        });
    });
}

function setupOpenCloseAll() {
    const openAll = document.getElementById("open-all");
    const closeAll = document.getElementById("close-all");

    if (openAll) {
        openAll.addEventListener("click", () => {
            document.querySelectorAll(".api-card").forEach((card) => {
                if (!card.hidden) {
                    card.open = true;
                }
            });
        });
    }

    if (closeAll) {
        closeAll.addEventListener("click", () => {
            document.querySelectorAll(".api-card").forEach((card) => {
                card.open = false;
            });
        });
    }
}

function setupSectionObserver() {
    const links = Array.from(document.querySelectorAll(".sidebar-nav a"));
    const targets = links
        .map((link) => document.querySelector(link.getAttribute("href")))
        .filter(Boolean);

    const observer = new IntersectionObserver(
        (entries) => {
            entries.forEach((entry) => {
                if (!entry.isIntersecting) {
                    return;
                }

                links.forEach((link) => {
                    const active = link.getAttribute("href") === `#${entry.target.id}`;
                    link.classList.toggle("active", active);
                });
            });
        },
        { rootMargin: "-35% 0px -55% 0px", threshold: 0.01 }
    );

    targets.forEach((section) => observer.observe(section));
}

function setupNavToggle() {
    const toggle = document.getElementById("nav-toggle");
    const sidebar = document.getElementById("sidebar");

    if (!toggle || !sidebar) {
        return;
    }

    toggle.addEventListener("click", () => {
        const open = sidebar.classList.toggle("open");
        toggle.setAttribute("aria-expanded", open ? "true" : "false");
    });

    document.querySelectorAll(".sidebar-nav a").forEach((link) => {
        link.addEventListener("click", () => {
            sidebar.classList.remove("open");
            toggle.setAttribute("aria-expanded", "false");
        });
    });
}

renderReference();
setupExampleTabs();
setupApiSearch();
setupOpenCloseAll();
setupSectionObserver();
setupNavToggle();
