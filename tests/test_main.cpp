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

int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;
    
    test_robust_qr_decomposition();
    test_eigenvalues_convergence();
    test_flattened_vector_invariants();
    test_pivot_selection_and_rank();
    test_precision_row_scaling();
    
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
#include <cassert>
#include <cmath>
#include "elda/matrix.hpp"

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

    // Dimensions check for Reduced QR: Q is 2x3, R is 3x3
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

    // 3. Complex spectrum matrix (90-degree rotation)
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

    // Verify safe accessor indexing mapped correctly
    assert(A(0, 0) == 1.0);
    assert(A(0, 1) == 2.0);
    assert(A(1, 0) == 3.0);
    assert(A(1, 1) == 4.0);

    // Verify cache locality arithmetic multiplication performance layout
    linalg::matrix B = A * A;
    // [1 2; 3 4] * [1 2; 3 4] = [7 10; 15 22]
    assert(is_close(B(0, 0), 7.0));
    assert(is_close(B(0, 1), 10.0));
    assert(is_close(B(1, 0), 15.0));
    assert(is_close(B(1, 1), 22.0));

    std::cout << "Flat vector layout invariants verified successfully!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;
    test_robust_qr_decomposition();
    test_eigenvalues_convergence();
    test_flattened_vector_invariants();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
#include <cassert>
#include <cmath>
#include <vector>
#include <stdexcept>
#include "elda/matrix.hpp"
#include <sstream>
#include <elda/linalg.hpp>

// Utility helper to assert floating point values within a safe numerical tolerance
bool is_close(double a, double b, double epsilon = 1e-4) {
    return std::abs(a - b) < epsilon;
}

// Helper lambda runner for verifying runtime exceptions
template <typename Func>
void expect_runtime_error(Func func) {
    try {
        func();
    } catch (const std::runtime_error&) {
        return;
    }
    assert(false && "Expected std::runtime_error was not thrown!");
}

void test_matrix_construction_and_identity() {
    std::cout << "Running construction and identity tests..." << std::endl;
    linalg::matrix m(2, 2);
    assert(m.get_rows() == 2);
    assert(m.get_cols() == 2);
    assert(is_close(m(0, 0), 0.0));

    linalg::matrix i2 = linalg::identity(2);
    assert(is_close(i2(0, 0), 1.0));
    assert(is_close(i2(0, 1), 0.0));
    assert(is_close(i2(1, 0), 0.0));
    assert(is_close(i2(1, 1), 1.0));

    std::cout << "Construction and identity tests passed!" << std::endl;
    std::cout << "Running construction tests..." << std::endl;

    linalg::matrix m(2, 2);
    assert(m.row == 2);
    assert(m.col == 2);
    assert(is_close(m.get_element(0, 0), 0.0));

    linalg::matrix filled(2, 2, 3.5);
    assert(is_close(filled.get_element(0, 0), 3.5));
    assert(is_close(filled.get_element(1, 1), 3.5));

    std::ostringstream out;
    out << filled;
    assert(out.str() == "3.5 3.5 \n3.5 3.5 \n");

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
    
    m1(0, 0) = 2.0; m1(0, 1) = 2.0;
    m1(1, 0) = 2.0; m1(1, 1) = 2.0;
    
    m2(0, 0) = 3.0; m2(0, 1) = 3.0;
    m2(1, 0) = 3.0; m2(1, 1) = 3.0;

    linalg::matrix sum = m1 + m2;
    assert(is_close(sum(0, 0), 5.0));
    assert(is_close(sum(1, 1), 5.0));

    linalg::matrix scaled = m1 * 2.0;
    assert(is_close(scaled(0, 0), 4.0));
    std::cout << "Arithmetic tests passed!" << std::endl;
}

void test_matrix_transpose_and_determinant() {
    std::cout << "Running transpose & determinant tests..." << std::endl;
    linalg::matrix m(2, 2);
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;

    assert(is_close(m.det(), -2.0));
    linalg::matrix t = m.transpose();
    assert(is_close(t(0, 1), 3.0));
    assert(is_close(t(1, 0), 2.0));
    std::cout << "Transpose & determinant tests passed!" << std::endl;
}

void test_inverse_rejects_singular_pivots() {
    std::cout << "Running inverse singular-pivot tests..." << std::endl;

    linalg::matrix zero(2, 2);
    expect_runtime_error([&zero]() { zero.inverse(); });

    linalg::matrix duplicate_rows(2, 2);
    duplicate_rows(0, 0) = 2.0; duplicate_rows(0, 1) = 4.0;
    duplicate_rows(1, 0) = 2.0; duplicate_rows(1, 1) = 4.0;
    expect_runtime_error([&duplicate_rows]() { duplicate_rows.inverse(); });

    linalg::matrix needs_row_swap(2, 2);
    needs_row_swap(0, 0) = 0.0; needs_row_swap(0, 1) = 2.0;
    needs_row_swap(1, 0) = 1.0; needs_row_swap(1, 1) = 3.0;
    
    linalg::matrix inverse = needs_row_swap.inverse();
    assert(is_close(inverse(0, 0), -1.5));
    assert(is_close(inverse(0, 1), 1.0));
    assert(is_close(inverse(1, 0), 0.5));
    assert(is_close(inverse(1, 1), 0.0));

    std::cout << "Inverse singular-pivot tests passed!" << std::endl;
}

void test_matpow_validation_and_results() {
    std::cout << "Running matpow tests..." << std::endl;
    linalg::matrix base(2, 2);
    base(0, 0) = 2.0; base(0, 1) = 1.0;
    base(1, 0) = 1.0; base(1, 1) = 2.0;

    linalg::matrix zero_power = linalg::matpow(base, 0);
    assert(is_close(zero_power(0, 0), 1.0));
    assert(is_close(zero_power(1, 1), 1.0));

    linalg::matrix squared = linalg::matpow(base, 2);
    assert(is_close(squared(0, 0), 5.0));
    assert(is_close(squared(0, 1), 4.0));

    std::cout << "Matpow tests passed!" << std::endl;
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
    linalg::matrix A_hat = Q * R;

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

    assert(Q_rect.get_rows() == 2 && Q_rect.get_cols() == 3);
    assert(R_rect.get_rows() == 3 && R_rect.get_cols() == 3);

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

    linalg::matrix D(2, 2);
    D(0, 0) = 5.0; D(0, 1) = 0.0;
    D(1, 0) = 0.0; D(1, 1) = -3.0;
    linalg::matrix evD = D.eigenvalues();
    assert(is_close(evD(0, 0), 5.0) || is_close(evD(0, 0), -3.0));

    linalg::matrix S(2, 2);
    S(0, 0) = 2.0; S(0, 1) = 1.0;
    S(1, 0) = 1.0; S(1, 1) = 2.0;
    linalg::matrix evS = S.eigenvalues();
    assert((is_close(evS(0, 0), 3.0) && is_close(evS(1, 0), 1.0)) || 
           (is_close(evS(0, 0), 1.0) && is_close(evS(1, 0), 3.0)));

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

    assert(A(0, 0) == 1.0);
    assert(A(0, 1) == 2.0);
    assert(A(1, 0) == 3.0);
    assert(A(1, 1) == 4.0);

    linalg::matrix B = A * A;
    assert(is_close(B(0, 0), 7.0));
    assert(is_close(B(0, 1), 10.0));
    assert(is_close(B(1, 0), 15.0));
    assert(is_close(B(1, 1), 22.0));

    std::cout << "Flat vector layout invariants verified successfully!" << std::endl;
}

void test_pivot_selection_and_rank() {
    std::cout << "Running Pivot Selection and Rank tests..." << std::endl;

    // 1. Matrix with leading zero column
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
    assert(T.rank() == 2);

    // 3. Wide matrix configuration (2x4)
    linalg::matrix W(2, 4);
    W(0, 0) = 1.0; W(0, 1) = 2.0; W(0, 2) = 3.0; W(0, 3) = 4.0;
    W(1, 0) = 2.0; W(1, 1) = 4.0; W(1, 2) = 6.0; W(1, 3) = 8.0;
    assert(W.rank() == 1);

    std::cout << "Pivot Selection and Rank tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;
    
    test_matrix_construction_and_identity();
    test_matrix_arithmetic();
    test_matrix_transpose_and_determinant();
    test_inverse_rejects_singular_pivots();
    test_matpow_validation_and_results();
    test_robust_qr_decomposition();
    test_eigenvalues_convergence();
    test_flattened_vector_invariants();
    test_pivot_selection_and_rank();
    
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
