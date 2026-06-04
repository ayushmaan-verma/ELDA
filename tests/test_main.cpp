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
    A.arr[0][0] = 1.0; A.arr[0][1] = 1.0; A.arr[0][2] = 2.0;
    A.arr[1][0] = 2.0; A.arr[1][1] = 2.0; A.arr[1][2] = 5.0;
    A.arr[2][0] = 3.0; A.arr[2][1] = 3.0; A.arr[2][2] = 8.0;

    linalg::matrix Q = A.qr_decomp_q();
    linalg::matrix R = A.qr_decomp_r();

    // Reconstruct A_hat = Q * R
    linalg::matrix A_hat = Q * R;

    // Verify reconstruction satisfies A == Q * R within tolerance
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            assert(is_close(A.arr[i][j], A_hat.arr[i][j]));
        }
    }

    // 2. Test case: Rectangular wide matrix (2 x 3)
    linalg::matrix B(2, 3);
    B.arr[0][0] = 1.0; B.arr[0][1] = 0.0; B.arr[0][2] = 2.0;
    B.arr[1][0] = 0.0; B.arr[1][1] = 1.0; B.arr[1][2] = 3.0;

    linalg::matrix Q_rect = B.qr_decomp_q();
    linalg::matrix R_rect = B.qr_decomp_r();

    // Dimensions check for Reduced QR: Q is 2x3, R is 3x3
    assert(Q_rect.row == 2 && Q_rect.col == 3);
    assert(R_rect.row == 3 && R_rect.col == 3);

    linalg::matrix B_hat = Q_rect * R_rect;
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 3; j++) {
            assert(is_close(B.arr[i][j], B_hat.arr[i][j]));
        }
    }

    std::cout << "Robust QR Decomposition tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;
    test_robust_qr_decomposition();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}