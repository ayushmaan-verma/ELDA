#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "elda/matrix.hpp"

bool is_close(double a, double b, double epsilon = 1e-4) {
    return std::abs(a - b) < epsilon;
}

void test_robust_qr_decomposition() {
    std::cout << "Running Robust QR Decomposition tests..." << std::endl;

    linalg::matrix A(3, 3);
    A(0, 0) = 1.0;
    A(0, 1) = 1.0;
    A(0, 2) = 2.0;
    A(1, 0) = 2.0;
    A(1, 1) = 2.0;
    A(1, 2) = 5.0;
    A(2, 0) = 3.0;
    A(2, 1) = 3.0;
    A(2, 2) = 8.0;

    linalg::matrix Q = A.qr_decomp_q();
    linalg::matrix R = A.qr_decomp_r();
    linalg::matrix A_hat = Q * R;

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            assert(is_close(A(i, j), A_hat(i, j)));
        }
    }

    linalg::matrix B(2, 3);
    B(0, 0) = 1.0;
    B(0, 1) = 0.0;
    B(0, 2) = 2.0;
    B(1, 0) = 0.0;
    B(1, 1) = 1.0;
    B(1, 2) = 3.0;

    linalg::matrix Q_rect = B.qr_decomp_q();
    linalg::matrix R_rect = B.qr_decomp_r();

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

    linalg::matrix D(2, 2);
    D(0, 0) = 5.0;
    D(0, 1) = 0.0;
    D(1, 0) = 0.0;
    D(1, 1) = -3.0;
    linalg::matrix evD = D.eigenvalues();
    assert(is_close(evD(0, 0), 5.0) || is_close(evD(0, 0), -3.0));

    linalg::matrix S(2, 2);
    S(0, 0) = 2.0;
    S(0, 1) = 1.0;
    S(1, 0) = 1.0;
    S(1, 1) = 2.0;
    linalg::matrix evS = S.eigenvalues();
    assert((is_close(evS(0, 0), 3.0) && is_close(evS(1, 0), 1.0)) ||
           (is_close(evS(0, 0), 1.0) && is_close(evS(1, 0), 3.0)));

    linalg::matrix C(2, 2);
    C(0, 0) = 0.0;
    C(0, 1) = -1.0;
    C(1, 0) = 1.0;
    C(1, 1) = 0.0;

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
    A(0, 0) = 1.0;
    A(0, 1) = 2.0;
    A(1, 0) = 3.0;
    A(1, 1) = 4.0;

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
    Z(0, 0) = 0.0;
    Z(0, 1) = 1.0;
    Z(1, 0) = 0.0;
    Z(1, 1) = 2.0;
    assert(Z.rank() == 1);

    linalg::matrix T(4, 2);
    T(0, 0) = 1.0;
    T(0, 1) = 2.0;
    T(1, 0) = 2.0;
    T(1, 1) = 4.0;
    T(2, 0) = 3.0;
    T(2, 1) = 6.0;
    T(3, 0) = 1.0;
    T(3, 1) = -1.0;
    assert(T.rank() == 2);

    linalg::matrix W(2, 4);
    W(0, 0) = 1.0;
    W(0, 1) = 2.0;
    W(0, 2) = 3.0;
    W(0, 3) = 4.0;
    W(1, 0) = 2.0;
    W(1, 1) = 4.0;
    W(1, 2) = 6.0;
    W(1, 3) = 8.0;
    assert(W.rank() == 1);

    std::cout << "Pivot Selection and Rank tests passed!" << std::endl;
}

void test_precision_row_scaling() {
    std::cout << "Running precision row scaling regression tests..." << std::endl;

    linalg::matrix A(2, 2);
    A(0, 0) = 1e-7;
    A(0, 1) = 1.0;
    A(1, 0) = 1.0;
    A(1, 1) = 1.0;

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
    A(0, 0) = 1.0;
    A(0, 1) = 0.0;
    A(0, 2) = 0.0;
    A(0, 3) = 0.0;
    A(1, 0) = 0.0;
    A(1, 1) = 0.0;
    A(1, 2) = 0.0;
    A(1, 3) = 0.0;

    linalg::matrix B(2, 4);
    B(0, 0) = 0.0;
    B(0, 1) = 1.0;
    B(0, 2) = 0.0;
    B(0, 3) = 0.0;
    B(1, 0) = 0.0;
    B(1, 1) = 0.0;
    B(1, 2) = 0.0;
    B(1, 3) = 0.0;

    assert(is_close(linalg::inner_product(A, B), 0.0));
    assert(is_close(linalg::angle(A, B), linalg::PI / 2.0));

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

void test_solve_unique_with_row_swaps() {
    std::cout << "Running unique solution requiring row swaps tests..." << std::endl;

    linalg::matrix A(2, 3);
    A(0, 0) = 0.0;
    A(0, 1) = 1.0;
    A(0, 2) = 2.0;
    A(1, 0) = 1.0;
    A(1, 1) = 0.0;
    A(1, 2) = 3.0;

    linalg::matrix sol = A.solve();
    assert(sol.get_rows() == 2 && sol.get_cols() == 1);
    assert(is_close(sol(0, 0), 3.0));
    assert(is_close(sol(1, 0), 2.0));

    std::cout << "Unique solution with row swaps tests passed!" << std::endl;
}

void test_solve_inconsistent_system() {
    std::cout << "Running inconsistent system tests..." << std::endl;

    linalg::matrix A(2, 3);
    A(0, 0) = 1.0;
    A(0, 1) = 2.0;
    A(0, 2) = 3.0;
    A(1, 0) = 2.0;
    A(1, 1) = 4.0;
    A(1, 2) = 7.0;

    bool caught_exception = false;
    try {
        A.solve();
    } catch (const std::runtime_error& e) {
        caught_exception = true;
        assert(std::string(e.what()) == "Inconsistent system: No solution exists.");
    }

    assert(caught_exception && "Should throw for inconsistent system.");
    std::cout << "Inconsistent system tests passed!" << std::endl;
}

void test_solve_underdetermined_system() {
    std::cout << "Running underdetermined system tests..." << std::endl;

    linalg::matrix A(2, 3);
    A(0, 0) = 1.0;
    A(0, 1) = 2.0;
    A(0, 2) = 3.0;
    A(1, 0) = 2.0;
    A(1, 1) = 4.0;
    A(1, 2) = 6.0;

    bool caught_exception = false;
    try {
        A.solve();
    } catch (const std::runtime_error& e) {
        caught_exception = true;
        assert(std::string(e.what()) ==
               "Singular/Underdetermined system: Infinitely many solutions exist.");
    }

    assert(caught_exception && "Should throw for underdetermined system.");
    std::cout << "Underdetermined system tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;

    test_robust_qr_decomposition();
    test_eigenvalues_convergence();
    test_flattened_vector_invariants();
    test_pivot_selection_and_rank();
    test_precision_row_scaling();
    test_rectangular_geometric_utilities();
    test_solve_unique_with_row_swaps();
    test_solve_inconsistent_system();
    test_solve_underdetermined_system();

    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
