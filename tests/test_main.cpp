#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <sstream>
#include "elda/matrix.hpp"

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

    linalg::matrix Z(2, 2);
    Z(0, 0) = 0.0; Z(0, 1) = 1.0;
    Z(1, 0) = 0.0; Z(1, 1) = 2.0;
    assert(Z.rank() == 1);

    linalg::matrix T(4, 2);
    T(0, 0) = 1.0; T(0, 1) = 2.0;
    T(1, 0) = 2.0; T(1, 1) = 4.0;
    T(2, 0) = 3.0; T(2, 1) = 6.0;
    T(3, 0) = 1.0; T(3, 1) = -1.0;
    assert(T.rank() == 2);

    linalg::matrix W(2, 4);
    W(0, 0) = 1.0; W(0, 1) = 2.0; W(0, 2) = 3.0; W(0, 3) = 4.0;
    W(1, 0) = 2.0; W(1, 1) = 4.0; W(1, 2) = 6.0; W(1, 3) = 8.0;
    assert(W.rank() == 1);

    std::cout << "Pivot Selection and Rank tests passed!" << std::endl;
}

void test_precision_row_scaling() {
    std::cout << "Running precision row scaling regression tests..." << std::endl;

    linalg::matrix A(2, 2);
    A(0, 0) = 1e-7; A(0, 1) = 1.0;
    A(1, 0) = 1.0;  A(1, 1) = 1.0;

    linalg::matrix invA = A.inverse();
    linalg::matrix I_check = A * invA;
    
    assert(is_close(I_check(0, 0), 1.0, 1e-9));
    assert(is_close(I_check(1, 1), 1.0, 1e-9));
    assert(is_close(I_check(0, 1), 0.0, 1e-9));
    assert(is_close(I_check(1, 0), 0.0, 1e-9));

    std::cout << "Precision row scaling tests passed!" << std::endl;
}

void test_rectangular_geometric_utilities() {
    std::cout << "Running rectangular geometry inner product and angle tests..." << std::endl;

    linalg::matrix A(2, 4);
    A(0, 0) = 1.0; A(0, 1) = 0.0; A(0, 2) = 0.0; A(0, 3) = 0.0;
    A(1, 0) = 0.0; A(1, 1) = 0.0; A(1, 2) = 0.0; A(1, 3) = 0.0;

    linalg::matrix B(2, 4);
    B(0, 0) = 0.0; B(0, 1) = 1.0; B(0, 2) = 0.0; B(0, 3) = 0.0;
    B(1, 0) = 0.0; B(1, 1) = 0.0; B(1, 2) = 0.0; B(1, 3) = 0.0;

    double dot_prod = linalg::inner_product(A, B);
    assert(is_close(dot_prod, 0.0));

    double theta = linalg::angle(A, B);
    assert(is_close(theta, linalg::PI / 2.0));

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

void test_to_string_formatting() {
    std::cout << "Running to_string formatting tests..." << std::endl;

    linalg::matrix A(2, 2);
    A(0, 0) = 1.0; A(0, 1) = 2.0;
    A(1, 0) = 3.0; A(1, 1) = 4.0;

    std::string out = A.to_string();
    std::string expected = "1 2\n3 4\n";
    
    assert(out == expected && "to_string output did not match expected formatting!");

    std::cout << "to_string formatting tests passed!" << std::endl;
}

void test_solve_unique_with_row_swaps() {
    std::cout << "Running unique solution requiring row swaps tests..." << std::endl;
    linalg::matrix A(2, 3);
    A(0, 0) = 0.0; A(0, 1) = 1.0; A(0, 2) = 2.0;
    A(1, 0) = 1.0; A(1, 1) = 0.0; A(1, 2) = 3.0;

    linalg::matrix sol = A.solve();
    assert(sol.get_rows() == 2 && sol.get_cols() == 1);
    assert(is_close(sol(0, 0), 3.0));
    assert(is_close(sol(1, 0), 2.0));
    std::cout << "Unique solution with row swaps tests passed!" << std::endl;
}

void test_solve_inconsistent_system() {
    std::cout << "Running inconsistent system tests..." << std::endl;
    linalg::matrix A(2, 3);
    A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 3.0;
    A(1, 0) = 2.0; A(1, 1) = 4.0; A(1, 2) = 7.0;

    bool caught_exception = false;
    try {
        A.solve();
    } catch (const std::runtime_error& e) {
        caught_exception = true;
        std::string msg = e.what();
        assert(msg == "Inconsistent system: No solution exists.");
    }
    assert(caught_exception && "Should throw for inconsistent system.");
    std::cout << "Inconsistent system tests passed!" << std::endl;
}

void test_solve_underdetermined_system() {
    std::cout << "Running underdetermined system tests..." << std::endl;
    linalg::matrix A(2, 3);
    A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 3.0;
    A(1, 0) = 2.0; A(1, 1) = 4.0; A(1, 2) = 6.0;

    bool caught_exception = false;
    try {
        A.solve();
    } catch (const std::runtime_error& e) {
        caught_exception = true;
        std::string msg = e.what();
        assert(msg == "Singular/Underdetermined system: Infinitely many solutions exist.");
    }
    assert(caught_exception && "Should throw for underdetermined system.");
    std::cout << "Underdetermined system tests passed!" << std::endl;
}

void test_fpg_behavior() {
    std::cout << "Running FPG behavior tests..." << std::endl;

    // Create a matrix with sub-EPS entries
    linalg::matrix A(2, 2);
    A(0, 0) = 1e-10; A(0, 1) = 0.0;
    A(1, 0) = 0.0;  A(1, 1) = 1e-10;

    linalg::matrix B(2, 2);
    B(0, 0) = 2e-10; B(0, 1) = 0.0;
    B(1, 0) = 0.0;  B(1, 1) = 2e-10;

    // Test operator+
    linalg::matrix C = A + B;
    // Values should be exactly 3e-10, NOT cleaned up to 0.0
    assert(C(0, 0) == 3e-10);
    assert(C(1, 1) == 3e-10);

    // Test operator-
    linalg::matrix D = B - A;
    // Values should be exactly 1e-10, NOT cleaned up to 0.0
    assert(D(0, 0) == 1e-10);
    assert(D(1, 1) == 1e-10);

    // Test operator=
    linalg::matrix E(2, 2);
    E = C;
    assert(E(0, 0) == 3e-10);
    assert(E(1, 1) == 3e-10);

    // Explicit call to fpg() should cleanup the values
    linalg::fpg(E);
    assert(E(0, 0) == 0.0);
    assert(E(1, 1) == 0.0);

    std::cout << "FPG behavior tests passed!" << std::endl;
}

void test_matrix_dimension_validation() {
    std::cout << "Running matrix dimension validation tests..." << std::endl;

    // Test negative dimensions in constructor
    expect_runtime_error([]() { linalg::matrix m(-1, 3); });
    expect_runtime_error([]() { linalg::matrix m(3, -1); });
    expect_runtime_error([]() { linalg::matrix m(-5, -5); });
    expect_runtime_error([]() { linalg::matrix m(-1, 3, 2.0); });
    expect_runtime_error([]() { linalg::matrix m(3, -1, 2.0); });

    // Test negative identity size
    expect_runtime_error([]() { linalg::identity(-1); });
    expect_runtime_error([]() { linalg::identity(-10); });

    // Test zero-sized matrices
    linalg::matrix m0(0, 0);
    assert(m0.get_rows() == 0);
    assert(m0.get_cols() == 0);
    assert(m0.arr.size() == 0);

    linalg::matrix m05(0, 5);
    assert(m05.get_rows() == 0);
    assert(m05.get_cols() == 5);
    assert(m05.arr.size() == 0);

    linalg::matrix m50(5, 0);
    assert(m50.get_rows() == 5);
    assert(m50.get_cols() == 0);
    assert(m50.arr.size() == 0);

    std::cout << "Matrix dimension validation tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;
    
    test_matrix_construction_and_identity();
    test_matrix_dimension_validation();
    test_matrix_arithmetic();
    test_matrix_transpose_and_determinant();
    test_inverse_rejects_singular_pivots();
    test_matpow_validation_and_results();
    test_robust_qr_decomposition();
    test_eigenvalues_convergence();
    test_flattened_vector_invariants();
    test_pivot_selection_and_rank();
    test_precision_row_scaling();
    test_rectangular_geometric_utilities();
    test_to_string_formatting();
    test_fpg_behavior();
    
    test_solve_unique_with_row_swaps();
    test_solve_inconsistent_system();
    test_solve_underdetermined_system();
    
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}