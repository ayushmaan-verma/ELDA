#include <iostream>
#include <cassert>
#include <cmath>
#include "elda/matrix.hpp"

// Utility helper to assert floating point values within a safe numerical tolerance
bool is_close(double a, double b, double epsilon = 1e-4) {
    return std::abs(a - b) < epsilon;
}

void test_robust_qr_decomposition() {
    std::cout << "Running Robust QR Decomposition tests..." << std::endl;

    // 1. Test case: Rank-deficient matrix (Columns 1 and 2 are identical)
    linalg::matrix A(3, 3);
    A(0, 0) = 1.0; A(0, 1) = 1.0; A(0, 2) = 2.0;
    A(1, 0) = 2.0; A(1, 1) = 2.0; A(1, 2) = 5.0;
    A(2, 0) = 3.0; A(2, 1) = 3.0; A(2, 2) = 8.0;

    linalg::matrix Q = A.qr_decomp_q();
    linalg::matrix R = A.qr_decomp_r();

    // Reconstruct A_hat = Q * R
    linalg::matrix A_hat = Q * R;

    // Verify reconstruction satisfies A == Q * R within tolerance
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            assert(is_close(A(i, j), A_hat(i, j)));
        }
    }

    // 2. Test case: Rectangular wide matrix (2 x 3)
    linalg::matrix B(2, 3);
    B(0, 0) = 1.0; B(0, 1) = 0.0; B(0, 2) = 2.0;
    B(1, 0) = 0.0; B(1, 1) = 1.0; B(1, 2) = 3.0;

    linalg::matrix Q_rect = B.qr_decomp_q();
    linalg::matrix R_rect = B.qr_decomp_r();

    // Dimensions check for Reduced/Thin QR: Q is 2x3, R is 3x3
    assert(Q_rect.row == 2 && Q_rect.col == 3);
    assert(R_rect.row == 3 && R_rect.col == 3);

    linalg::matrix B_hat = Q_rect * R_rect;
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 3; j++) {
            assert(is_close(B(i, j), B_hat(i, j)));
        }
    }

    std::cout << "Robust QR Decomposition tests passed!" << std::endl;
}

void test_eigenvalues_convergence() {
    std::cout << "Running eigenvalue convergence tests..." << std::endl;

    // 1. Diagonal matrix (Immediate convergence)
    linalg::matrix D(2, 2);
    D(0, 0) = 5.0; D(0, 1) = 0.0;
    D(1, 0) = 0.0; D(1, 1) = -3.0;
    linalg::matrix evD = D.eigenvalues();
    assert(is_close(evD(0, 0), 5.0) || is_close(evD(0, 0), -3.0));

    // 2. Symmetric matrix (Guaranteed real eigenvalues)
    linalg::matrix S(2, 2);
    S(0, 0) = 2.0; S(0, 1) = 1.0;
    S(1, 0) = 1.0; S(1, 1) = 2.0;
    linalg::matrix evS = S.eigenvalues();
    assert((is_close(evS(0, 0), 3.0) && is_close(evS(1, 0), 1.0)) || 
           (is_close(evS(0, 0), 1.0) && is_close(evS(1, 0), 3.0)));

    // 3. Complex spectrum matrix (90-degree rotation matrix)
    linalg::matrix C(2, 2);
    C(0, 0) = 0.0;  C(0, 1) = -1.0;
    C(1, 0) = 1.0;  C(1, 1) = 0.0;
    
    bool caught_non_convergence = false;
    try {
        C.eigenvalues();
    } catch (const std::runtime_error&) {
        caught_non_convergence = true;
    }
    assert(caught_non_convergence && "Should throw on non-converging complex spectra.");

    std::cout << "Eigenvalue convergence tests passed!" << std::endl;
}

void test_flattened_vector_invariants() {
    std::cout << "Running flat vector representation tests..." << std::endl;

    linalg::matrix A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 2.0;
    A(1, 0) = 3.0; A(1, 1) = 4.0;

    // Verify safe accessor indexing mapped correctly on the contiguous layout
    assert(A(0, 0) == 1.0);
    assert(A(0, 1) == 2.0);
    assert(A(1, 0) == 3.0);
    assert(A(1, 1) == 4.0);

    linalg::matrix B = A * A;
    // Layout multiplication integrity verify: [1 2; 3 4]^2 = [7 10; 15 22]
    assert(is_close(B(0, 0), 7.0));
    assert(is_close(B(0, 1), 10.0));
    assert(is_close(B(1, 0), 15.0));
    assert(is_close(B(1, 1), 22.0));

    std::cout << "Flat vector layout invariants verified successfully!" << std::endl;
}

void test_pivot_selection_and_rank() {
    std::cout << "Running Pivot Selection and Rank tests..." << std::endl;

    // 1. Matrix with leading zero column (Rank deficient)
    // [0 1]
    // [0 2] -> rank should correctly evaluate to 1
    linalg::matrix Z(2, 2);
    Z(0, 0) = 0.0; Z(0, 1) = 1.0;
    Z(1, 0) = 0.0; Z(1, 1) = 2.0;
    assert(Z.rank() == 1);

    // 2. Tall matrix configuration (4x2)
    linalg::matrix T(4, 2);
    T(0, 0) = 1.0; T(0, 1) = 2.0;
    T(1, 0) = 2.0; T(1, 1) = 4.0;
    T(2, 0) = 3.0; T(2, 1) = 6.0;
    T(3, 0) = 1.0; T(3, 1) = -1.0;
    // Col 1 and Col 2 are linearly independent due to row 4 -> rank 2
    assert(T.rank() == 2);

    // 3. Wide matrix configuration (2x4)
    linalg::matrix W(2, 4);
    W(0, 0) = 1.0; W(0, 1) = 2.0; W(0, 2) = 3.0; W(0, 3) = 4.0;
    W(1, 0) = 2.0; W(1, 1) = 4.0; W(1, 2) = 6.0; W(1, 3) = 8.0;
    // Row 2 is a direct scalar multiple of Row 1 -> rank 1
    assert(W.rank() == 1);

    std::cout << "Pivot Selection and Rank tests passed!" << std::endl;
}
void test_precision_row_scaling() {
    std::cout << "Running precision row scaling regression tests..." << std::endl;

    // Ill-conditioned matrix with an extremely small pivot (1e-7)
    // If we use `1.0 / 1e-7`, IEEE 754 truncation artifacts can ruin the exactness.
    // Direct division preserves the stability perfectly.
    linalg::matrix A(2, 2);
    A(0,0) = 1e-7; A(0,1) = 1.0;
    A(1,0) = 1.0;  A(1,1) = 1.0;

    linalg::matrix invA = A.inverse();
    
    // Check structural identity constraint: A * invA == Identity Matrix
    linalg::matrix I_check = A * invA;
    
    assert(is_close(I_check(0,0), 1.0, 1e-9));
    assert(is_close(I_check(1,1), 1.0, 1e-9));
    assert(is_close(I_check(0,1), 0.0, 1e-9));
    assert(is_close(I_check(1,0), 0.0, 1e-9));

    std::cout << "Precision row scaling tests passed!" << std::endl;
}
void test_rectangular_geometric_utilities() {
    std::cout << "Running rectangular geometry inner product and angle tests..." << std::endl;

    // Create matching rectangular wide (2x4) matrices
    linalg::matrix A(2, 4);
    A(0,0) = 1.0; A(0,1) = 0.0; A(0,2) = 0.0; A(0,3) = 0.0;
    A(1,0) = 0.0; A(1,1) = 0.0; A(1,2) = 0.0; A(1,3) = 0.0;

    linalg::matrix B(2, 4);
    B(0,0) = 0.0; B(0,1) = 1.0; B(0,2) = 0.0; B(0,3) = 0.0;
    B(1,0) = 0.0; B(1,1) = 0.0; B(1,2) = 0.0; B(1,3) = 0.0;

    // These matrices are completely orthogonal -> inner product must be 0
    double dot_prod = linalg::inner_product(A, B);
    assert(is_close(dot_prod, 0.0));

    // Orthogonal spaces have an angle deviation of exactly Pi / 2 (90 degrees)
    double theta = linalg::angle(A, B);
    assert(is_close(theta, linalg::PI / 2.0));

    // Ensure it catches structural size mismatches defensively
    linalg::matrix mismatched(3, 3);
    bool caught_size_error = false;
    try {
        linalg::inner_product(A, mismatched);
    } catch (const std::runtime_error&) {
        caught_size_error = true;
    }
    assert(caught_size_error && "Should throw error if space layouts do not match perfectly.");

    std::cout << "Rectangular geometric utility tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;
    
    test_robust_qr_decomposition();
    test_eigenvalues_convergence();
    test_flattened_vector_invariants();
    test_pivot_selection_and_rank();
    test_precision_row_scaling();
    test_rectangular_geometric_utilities();
    
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <elda/linalg.hpp>

bool is_close(double a, double b, double epsilon = 1e-5) {
    return std::abs(a - b) < epsilon;
}

template <typename Func>
void expect_runtime_error(Func func) {
    try {
        func();
    } catch (const std::runtime_error&) {
        return;
    }
    assert(false && "Expected std::runtime_error");
}

void test_matrix_construction_and_identity() {
    std::cout << "Running construction tests..." << std::endl;
    linalg::matrix m(2, 2);
    assert(m.row == 2);
    assert(m.col == 2);
    assert(is_close(m.get_element(0, 0), 0.0));

    linalg::matrix filled(2, 2, 3.5);
    assert(is_close(filled.get_element(0, 0), 3.5));
    assert(is_close(filled.get_element(1, 1), 3.5));

    linalg::matrix empty(0, 0);
    assert(empty.row == 0);
    assert(empty.col == 0);
    assert(empty.arr.empty());

    linalg::matrix i2 = linalg::identity(2);
    assert(is_close(i2.arr[0][0], 1.0));
    assert(is_close(i2.arr[0][1], 0.0));
    assert(is_close(i2.arr[1][0], 0.0));
    assert(is_close(i2.arr[1][1], 1.0));

    expect_runtime_error([]() { linalg::matrix(-1, 2); });
    expect_runtime_error([]() { linalg::matrix(2, -1); });
    expect_runtime_error([]() { linalg::identity(-1); });
    std::cout << "Construction tests passed!" << std::endl;
}

void test_matrix_arithmetic() {
    std::cout << "Running arithmetic tests..." << std::endl;
    linalg::matrix m1(2, 2);
    linalg::matrix m2(2, 2);
    m1.arr = {{2.0, 2.0}, {2.0, 2.0}};
    m2.arr = {{3.0, 3.0}, {3.0, 3.0}};
    linalg::matrix sum = m1 + m2;
    assert(is_close(sum.get_element(0, 0), 5.0));
    assert(is_close(sum.get_element(1, 1), 5.0));
    linalg::matrix scaled = m1 * 2.0;
    assert(is_close(scaled.get_element(0, 0), 4.0));
    std::cout << "Arithmetic tests passed!" << std::endl;
}

void test_matrix_transpose_and_determinant() {
    std::cout << "Running transpose & determinant tests..." << std::endl;
    linalg::matrix m(2, 2);
    m.arr = {{1.0, 2.0}, {3.0, 4.0}};
    assert(is_close(m.det(), -2.0));
    linalg::matrix t = m.transpose();
    assert(is_close(t.get_element(0, 1), 3.0));
    assert(is_close(t.get_element(1, 0), 2.0));
    std::cout << "Transpose & determinant tests passed!" << std::endl;
}

void test_inverse_rejects_singular_pivots() {
    std::cout << "Running inverse singular-pivot tests..." << std::endl;

    linalg::matrix zero(2, 2);
    expect_runtime_error([&zero]() { zero.inverse(); });

    linalg::matrix duplicate_rows(2, 2);
    duplicate_rows.arr = {{2.0, 4.0}, {2.0, 4.0}};
    expect_runtime_error([&duplicate_rows]() { duplicate_rows.inverse(); });

    linalg::matrix near_singular(2, 2);
    near_singular.arr = {{1.0, 1.0}, {1.0, 1.0 + (linalg::EPS / 2.0)}};
    expect_runtime_error([&near_singular]() { near_singular.inverse(); });

    linalg::matrix needs_row_swap(2, 2);
    needs_row_swap.arr = {{0.0, 2.0}, {1.0, 3.0}};
    linalg::matrix inverse = needs_row_swap.inverse();
    assert(is_close(inverse.get_element(0, 0), -1.5));
    assert(is_close(inverse.get_element(0, 1), 1.0));
    assert(is_close(inverse.get_element(1, 0), 0.5));
    assert(is_close(inverse.get_element(1, 1), 0.0));

    std::cout << "Inverse singular-pivot tests passed!" << std::endl;
}

void test_transforms_and_vectors() {
    std::cout << "Running transforms & vector utilities tests..." << std::endl;
    linalg::matrix trans = linalg::shift(5.0, 10.0);
    assert(is_close(trans.get_element(0, 2), 5.0));
    assert(is_close(trans.get_element(1, 2), 10.0));
    linalg::matrix v3 = linalg::vec3(4.5, 5.5, 6.5);
    assert(v3.row == 3);
    assert(v3.col == 1);
    assert(is_close(v3.get_element(0, 0), 4.5));
    assert(is_close(v3.get_element(2, 0), 6.5));
    std::cout << "Transforms & vector utilities tests passed!" << std::endl;
}

void test_vector_shape_validation() {
    std::cout << "Running check_lin_comb validation tests..." << std::endl;
    linalg::matrix v1(3, 1, 1.0);
    linalg::matrix v2(3, 1, 2.0);
    linalg::matrix target(3, 1, 5.0);

    linalg::matrix bad_row_vec(2, 1, 1.0);
    bool caught_row_mismatch = false;
    try {
        std::vector<linalg::matrix> vec_set = {v1, bad_row_vec};
        linalg::check_lin_comb(target, vec_set);
    } catch (const std::runtime_error& e) {
        caught_row_mismatch = true;
    }
    assert(caught_row_mismatch && "Should throw exception for mismatched row sizes");

    linalg::matrix wide_matrix(3, 2, 1.0);
    bool caught_non_column = false;
    try {
        std::vector<linalg::matrix> vec_set = {v1, wide_matrix};
        linalg::check_lin_comb(target, vec_set);
    } catch (const std::runtime_error& e) {
        caught_non_column = true;
    }
    assert(caught_non_column && "Should throw exception if any input is not a column vector");

    std::cout << "check_lin_comb validation tests passed!" << std::endl;
}

void test_matpow_validation_and_results() {
    std::cout << "Running matpow tests..." << std::endl;
    linalg::matrix base(2, 2);
    base.arr = {{2.0, 1.0}, {1.0, 2.0}};

    linalg::matrix zero_power = linalg::matpow(base, 0);
    assert(is_close(zero_power.get_element(0, 0), 1.0));
    assert(is_close(zero_power.get_element(0, 1), 0.0));
    assert(is_close(zero_power.get_element(1, 0), 0.0));
    assert(is_close(zero_power.get_element(1, 1), 1.0));

    linalg::matrix first_power = linalg::matpow(base, 1);
    assert(first_power == base);

    linalg::matrix squared = linalg::matpow(base, 2);
    assert(is_close(squared.get_element(0, 0), 5.0));
    assert(is_close(squared.get_element(0, 1), 4.0));
    assert(is_close(squared.get_element(1, 0), 4.0));
    assert(is_close(squared.get_element(1, 1), 5.0));

    bool rejected_negative_exponent = false;
    try {
        linalg::matpow(base, -1);
    } catch (const std::runtime_error&) {
        rejected_negative_exponent = true;
    }
    assert(rejected_negative_exponent);

    bool rejected_non_square_matrix = false;
    try {
        linalg::matrix non_square(2, 3);
        linalg::matpow(non_square, 0);
    } catch (const std::runtime_error&) {
        rejected_non_square_matrix = true;
    }
    assert(rejected_non_square_matrix);

    std::cout << "Matpow tests passed!" << std::endl;
}

bool matrix_is_finite(linalg::matrix m) {
    for (int i = 0; i < m.row; i++) {
        for (int j = 0; j < m.col; j++) {
            if (!std::isfinite(m.arr[i][j])) {
                return false;
            }
        }
    }
    return true;
}

void assert_reconstructs(linalg::matrix original) {
    linalg::matrix q = original.qr_decomp_q();
    linalg::matrix r = original.qr_decomp_r();
    linalg::matrix rebuilt = q * r;

    assert(q.row == original.row);
    assert(q.col == original.col);
    assert(r.row == original.col);
    assert(r.col == original.col);
    assert(matrix_is_finite(q));
    assert(matrix_is_finite(r));
    assert(matrix_is_finite(rebuilt));

    for (int i = 0; i < original.row; i++) {
        for (int j = 0; j < original.col; j++) {
            assert(is_close(original.arr[i][j], rebuilt.arr[i][j], 1e-4));
        }
    }
}

void test_gram_schmidt_bounds_and_dependent_columns() {
    std::cout << "Running Gram-Schmidt bounds and dependency tests..." << std::endl;

    linalg::matrix one_column(3, 1);
    one_column.arr = {{3.0}, {4.0}, {0.0}};
    linalg::matrix one_column_q = one_column.orthonormalize();
    assert(matrix_is_finite(one_column_q));
    assert(is_close(one_column_q.arr[0][0], 0.6));
    assert(is_close(one_column_q.arr[1][0], 0.8));
    assert_reconstructs(one_column);

    linalg::matrix rectangular(2, 3);
    rectangular.arr = {{1.0, 0.0, 2.0}, {0.0, 1.0, 3.0}};
    assert_reconstructs(rectangular);

    linalg::matrix zero_column(3, 2);
    zero_column.arr = {{0.0, 1.0}, {0.0, 2.0}, {0.0, 3.0}};
    linalg::matrix zero_column_q = zero_column.orthonormalize();
    assert(matrix_is_finite(zero_column_q));
    assert(is_close(zero_column_q.arr[0][0], 0.0));
    assert(is_close(zero_column_q.arr[1][0], 0.0));
    assert(is_close(zero_column_q.arr[2][0], 0.0));
    assert_reconstructs(zero_column);

    linalg::matrix dependent(3, 3);
    dependent.arr = {{1.0, 1.0, 2.0}, {2.0, 2.0, 5.0}, {3.0, 3.0, 8.0}};
    linalg::matrix dependent_q = dependent.orthonormalize();
    assert(matrix_is_finite(dependent_q));
    assert(is_close(dependent_q.arr[0][1], 0.0));
    assert(is_close(dependent_q.arr[1][1], 0.0));
    assert(is_close(dependent_q.arr[2][1], 0.0));
    assert_reconstructs(dependent);

    std::cout << "Gram-Schmidt bounds and dependency tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA UNIT TESTS ===" << std::endl;
    test_matrix_construction_and_identity();
    test_matrix_arithmetic();
    test_matrix_transpose_and_determinant();
    test_inverse_rejects_singular_pivots();
    test_transforms_and_vectors();
    test_vector_shape_validation();
    test_matpow_validation_and_results();
    test_gram_schmidt_bounds_and_dependent_columns();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
