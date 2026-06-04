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

void test_eigenvalues_convergence() {
    std::cout << "Running eigenvalue convergence tests..." << std::endl;

    // 1. Diagonal matrix (Immediate convergence)
    linalg::matrix D(2, 2);
    D.arr[0][0] = 5.0; D.arr[0][1] = 0.0;
    D.arr[1][0] = 0.0; D.arr[1][1] = -3.0;
    linalg::matrix evD = D.eigenvalues();
    assert(is_close(evD.arr[0][0], 5.0) || is_close(evD.arr[0][0], -3.0));

    // 2. Symmetric matrix (Guaranteed real eigenvalues)
    linalg::matrix S(2, 2);
    S.arr[0][0] = 2.0; S.arr[0][1] = 1.0;
    S.arr[1][0] = 1.0; S.arr[1][1] = 2.0;
    linalg::matrix evS = S.eigenvalues();
    assert((is_close(evS.arr[0][0], 3.0) && is_close(evS.arr[1][0], 1.0)) || 
           (is_close(evS.arr[0][0], 1.0) && is_close(evS.arr[1][0], 3.0)));

    // 3. Complex spectrum matrix (90-degree rotation)
    linalg::matrix C(2, 2);
    C.arr[0][0] = 0.0;  C.arr[0][1] = -1.0;
    C.arr[1][0] = 1.0;  C.arr[1][1] = 0.0;
    
    bool caught_non_convergence = false;
    try {
        C.eigenvalues();
    } catch (const std::runtime_error&) {
        caught_non_convergence = true;
    }
    assert(caught_non_convergence && "Should throw on non-converging complex spectra.");

    std::cout << "Eigenvalue convergence tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;
    test_robust_qr_decomposition();
    test_eigenvalues_convergence();
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

int main() {
    std::cout << "=== STARTING ELDA UNIT TESTS ===" << std::endl;
    test_matrix_construction_and_identity();
    test_matrix_arithmetic();
    test_matrix_transpose_and_determinant();
    test_transforms_and_vectors();
    test_vector_shape_validation();
    test_matpow_validation_and_results();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
