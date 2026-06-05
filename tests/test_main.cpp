#include <iostream>
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