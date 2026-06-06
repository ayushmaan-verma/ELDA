// tests/test_main.cpp
// ELDA core unit tests – single, authoritative entry point.
//
// New in this revision (issue #6 – tolerant floating-point comparisons):
//   test_approx_equal_helper          – scalar helper edge-cases
//   test_matrices_approx_equal        – matrix-level helper + shape mismatch
//   test_check_ortho_rotation         – rotation matrices via std::cos/sin
//   test_check_ortho_perturbed        – tiny perturbations that exact == rejects
//   test_triangular_checks_tolerance  – near-zero sub-/super-diagonal entries

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>
#include "elda/matrix.hpp"

// ---------------------------------------------------------------------------
// Scalar tolerance helper (independent of the library's approx_equal so that
// tests remain self-contained even if the library function is broken).
// ---------------------------------------------------------------------------
static bool is_close(double a, double b, double eps = 1e-6) {
    return std::abs(a - b) <= eps;
}

// Runs func(); asserts that it throws std::runtime_error.
template <typename Func>
static void expect_runtime_error(Func func) {
    try {
        func();
    } catch (const std::runtime_error&) {
        return;
    }
    assert(false && "Expected std::runtime_error was not thrown!");
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a 2x2 rotation matrix for angle theta (radians).
static linalg::matrix rotation2d(double theta) {
    linalg::matrix R(2, 2);
    R(0, 0) =  std::cos(theta); R(0, 1) = -std::sin(theta);
    R(1, 0) =  std::sin(theta); R(1, 1) =  std::cos(theta);
    return R;
}

// ---------------------------------------------------------------------------
// Issue #6 – new tests
// ---------------------------------------------------------------------------

void test_approx_equal_helper() {
    std::cout << "Running approx_equal scalar helper tests..." << std::endl;

    // Exact equality
    assert(linalg::approx_equal(1.0, 1.0));
    assert(linalg::approx_equal(0.0, 0.0));

    // Within default tolerance (1e-9)
    assert(linalg::approx_equal(1.0, 1.0 + 1e-10));
    assert(linalg::approx_equal(0.0, 1e-10));

    // Just at the boundary
    assert(linalg::approx_equal(0.0, linalg::APPROX_TOL));   // exactly on boundary → true
    assert(!linalg::approx_equal(0.0, linalg::APPROX_TOL * 2.0)); // outside → false

    // Negative-zero and positive-zero are equal
    assert(linalg::approx_equal(-0.0, 0.0));

    // Custom tolerance
    assert( linalg::approx_equal(1.0, 1.001, 1e-2));
    assert(!linalg::approx_equal(1.0, 1.001, 1e-4));

    std::cout << "approx_equal scalar helper tests passed!" << std::endl;
}

void test_matrices_approx_equal() {
    std::cout << "Running matrices_approx_equal tests..." << std::endl;

    linalg::matrix A(2, 2);
    A(0,0) = 1.0; A(0,1) = 2.0;
    A(1,0) = 3.0; A(1,1) = 4.0;

    linalg::matrix B(2, 2);
    B(0,0) = 1.0 + 1e-11; B(0,1) = 2.0 - 1e-11;
    B(1,0) = 3.0;         B(1,1) = 4.0 + 5e-12;

    // Within default APPROX_TOL
    assert( linalg::matrices_approx_equal(A, B));

    // Shape mismatch: must return false immediately (no access to element storage)
    linalg::matrix C(3, 3);
    assert(!linalg::matrices_approx_equal(A, C));

    // Values outside tolerance
    linalg::matrix D(2, 2);
    D(0,0) = 1.0 + 1e-7; D(0,1) = 2.0;
    D(1,0) = 3.0;        D(1,1) = 4.0;
    assert(!linalg::matrices_approx_equal(A, D, 1e-9));
    // But passes with looser tolerance
    assert( linalg::matrices_approx_equal(A, D, 1e-6));

    std::cout << "matrices_approx_equal tests passed!" << std::endl;
}

void test_check_ortho_rotation() {
    std::cout << "Running check_ortho rotation-matrix tests..." << std::endl;

    // 30-degree, 45-degree, and 90-degree 2×2 rotation matrices are all orthogonal.
    const double angles[] = {
        linalg::PI / 6.0,   // 30°
        linalg::PI / 4.0,   // 45°
        linalg::PI / 3.0,   // 60°
        linalg::PI / 2.0,   // 90°
        linalg::PI,         // 180°
    };

    for (double theta : angles) {
        linalg::matrix R = rotation2d(theta);
        // Must pass with default tolerance even though cos/sin introduce rounding.
        assert(linalg::check_ortho(R) &&
               "Rotation matrix should be detected as orthogonal");
    }

    // A non-orthogonal matrix must fail.
    linalg::matrix nonOrtho(2, 2);
    nonOrtho(0,0) = 2.0; nonOrtho(0,1) = 0.0;
    nonOrtho(1,0) = 0.0; nonOrtho(1,1) = 1.0;
    assert(!linalg::check_ortho(nonOrtho) &&
           "Non-orthogonal matrix must not pass check_ortho");

    // Singular matrix (no inverse) must return false, not throw.
    linalg::matrix singular(2, 2);
    singular(0,0) = 1.0; singular(0,1) = 2.0;
    singular(1,0) = 2.0; singular(1,1) = 4.0;
    assert(!linalg::check_ortho(singular));

    // Non-square matrix must return false.
    linalg::matrix rect(2, 3);
    assert(!linalg::check_ortho(rect));

    std::cout << "check_ortho rotation-matrix tests passed!" << std::endl;
}

void test_check_ortho_perturbed() {
    std::cout << "Running check_ortho tiny-perturbation tests..." << std::endl;

    // Build exact 45-degree rotation then add a perturbation smaller than
    // EPS (1e-6).  The result is still mathematically orthogonal to machine
    // precision and must pass with the default tolerance.
    linalg::matrix R = rotation2d(linalg::PI / 4.0);
    // Add perturbation ~ 1e-12 on one entry (much less than EPS = 1e-6).
    double saved = R(0, 0);
    R(0, 0) = saved + 1e-12;

    assert(linalg::check_ortho(R) &&
           "Rotation matrix with sub-tolerance perturbation should still pass");

    // A perturbation larger than EPS should fail with a tight explicit tol.
    linalg::matrix R2 = rotation2d(linalg::PI / 4.0);
    R2(0, 0) = R2(0, 0) + 1e-4;
    assert(!linalg::check_ortho(R2, 1e-9) &&
           "Large perturbation must fail strict tolerance check");

    std::cout << "check_ortho perturbation tests passed!" << std::endl;
}

void test_triangular_checks_tolerance() {
    std::cout << "Running triangular-check tolerance tests..." << std::endl;

    // Upper-triangular: sub-diagonal entries are tiny rounding residuals.
    linalg::matrix U(3, 3);
    U(0,0) = 1.0; U(0,1) = 2.0; U(0,2) = 3.0;
    U(1,0) = 5e-11; U(1,1) = 4.0; U(1,2) = 5.0; // tiny sub-diagonal residual
    U(2,0) = 0.0;   U(2,1) = 3e-12; U(2,2) = 6.0;
    assert(U.check_upper_tri() &&
           "Near-zero sub-diagonal should be accepted by tolerant check_upper_tri");

    // Genuinely non-upper-triangular entry must still fail.
    U(2, 0) = 0.5;
    assert(!U.check_upper_tri());

    // Lower-triangular: super-diagonal entries are tiny rounding residuals.
    linalg::matrix L(3, 3);
    L(0,0) = 1.0; L(0,1) = 4e-12; L(0,2) = 0.0; // tiny super-diagonal residual
    L(1,0) = 2.0; L(1,1) = 5.0;   L(1,2) = 7e-11;
    L(2,0) = 3.0; L(2,1) = 6.0;   L(2,2) = 9.0;
    assert(L.check_lower_tri() &&
           "Near-zero super-diagonal should be accepted by tolerant check_lower_tri");

    L(0, 2) = 0.5;
    assert(!L.check_lower_tri());

    std::cout << "Triangular-check tolerance tests passed!" << std::endl;
}

// ---------------------------------------------------------------------------
// Existing tests (kept from prior history, using operator() accessor)
// ---------------------------------------------------------------------------

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
    m1(0,0) = 2.0; m1(0,1) = 2.0;
    m1(1,0) = 2.0; m1(1,1) = 2.0;
    m2(0,0) = 3.0; m2(0,1) = 3.0;
    m2(1,0) = 3.0; m2(1,1) = 3.0;

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
    m(0,0) = 1.0; m(0,1) = 2.0;
    m(1,0) = 3.0; m(1,1) = 4.0;

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

    linalg::matrix dup(2, 2);
    dup(0,0) = 2.0; dup(0,1) = 4.0;
    dup(1,0) = 2.0; dup(1,1) = 4.0;
    expect_runtime_error([&dup]() { dup.inverse(); });

    linalg::matrix nrs(2, 2);
    nrs(0,0) = 0.0; nrs(0,1) = 2.0;
    nrs(1,0) = 1.0; nrs(1,1) = 3.0;
    linalg::matrix inv = nrs.inverse();
    assert(is_close(inv(0, 0), -1.5));
    assert(is_close(inv(0, 1),  1.0));
    assert(is_close(inv(1, 0),  0.5));
    assert(is_close(inv(1, 1),  0.0));

    std::cout << "Inverse singular-pivot tests passed!" << std::endl;
}

void test_matpow_validation_and_results() {
    std::cout << "Running matpow tests..." << std::endl;

    linalg::matrix base(2, 2);
    base(0,0) = 2.0; base(0,1) = 1.0;
    base(1,0) = 1.0; base(1,1) = 2.0;

    linalg::matrix zp = linalg::matpow(base, 0);
    assert(is_close(zp(0,0), 1.0));
    assert(is_close(zp(1,1), 1.0));

    linalg::matrix sq = linalg::matpow(base, 2);
    assert(is_close(sq(0,0), 5.0));
    assert(is_close(sq(0,1), 4.0));

    std::cout << "Matpow tests passed!" << std::endl;
}

void test_robust_qr_decomposition() {
    std::cout << "Running Robust QR Decomposition tests..." << std::endl;

    // 1. Rank-deficient 3x3
    linalg::matrix A(3, 3);
    A(0,0) = 1.0; A(0,1) = 1.0; A(0,2) = 2.0;
    A(1,0) = 2.0; A(1,1) = 2.0; A(1,2) = 5.0;
    A(2,0) = 3.0; A(2,1) = 3.0; A(2,2) = 8.0;

    linalg::matrix Q = A.qr_decomp_q();
    linalg::matrix R = A.qr_decomp_r();
    linalg::matrix Ah = Q * R;
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
            assert(is_close(A(i,j), Ah(i,j)));

    // 2. Rectangular 2x3
    linalg::matrix B(2, 3);
    B(0,0) = 1.0; B(0,1) = 0.0; B(0,2) = 2.0;
    B(1,0) = 0.0; B(1,1) = 1.0; B(1,2) = 3.0;

    linalg::matrix Qr = B.qr_decomp_q();
    linalg::matrix Rr = B.qr_decomp_r();
    assert(Qr.get_rows() == 2 && Qr.get_cols() == 3);
    assert(Rr.get_rows() == 3 && Rr.get_cols() == 3);

    linalg::matrix Bh = Qr * Rr;
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 3; j++)
            assert(is_close(B(i,j), Bh(i,j)));

    std::cout << "Robust QR Decomposition tests passed!" << std::endl;
}

void test_eigenvalues_convergence() {
    std::cout << "Running eigenvalue convergence tests..." << std::endl;

    linalg::matrix D(2, 2);
    D(0,0) = 5.0; D(0,1) = 0.0;
    D(1,0) = 0.0; D(1,1) = -3.0;
    linalg::matrix evD = D.eigenvalues();
    assert(is_close(evD(0,0), 5.0) || is_close(evD(0,0), -3.0));

    linalg::matrix S(2, 2);
    S(0,0) = 2.0; S(0,1) = 1.0;
    S(1,0) = 1.0; S(1,1) = 2.0;
    linalg::matrix evS = S.eigenvalues();
    assert((is_close(evS(0,0), 3.0) && is_close(evS(1,0), 1.0)) ||
           (is_close(evS(0,0), 1.0) && is_close(evS(1,0), 3.0)));

    linalg::matrix C(2, 2);
    C(0,0) = 0.0; C(0,1) = -1.0;
    C(1,0) = 1.0; C(1,1) =  0.0;
    bool caught = false;
    try { C.eigenvalues(); } catch (const std::runtime_error&) { caught = true; }
    assert(caught && "Should throw on non-converging complex spectra.");

    std::cout << "Eigenvalue convergence tests passed!" << std::endl;
}

void test_flattened_vector_invariants() {
    std::cout << "Running flat vector representation tests..." << std::endl;

    linalg::matrix A(2, 2);
    A(0,0) = 1.0; A(0,1) = 2.0;
    A(1,0) = 3.0; A(1,1) = 4.0;

    assert(A(0,0) == 1.0);
    assert(A(0,1) == 2.0);
    assert(A(1,0) == 3.0);
    assert(A(1,1) == 4.0);

    linalg::matrix B = A * A; // [7 10; 15 22]
    assert(is_close(B(0,0),  7.0));
    assert(is_close(B(0,1), 10.0));
    assert(is_close(B(1,0), 15.0));
    assert(is_close(B(1,1), 22.0));

    std::cout << "Flat vector layout invariants verified successfully!" << std::endl;
}

void test_pivot_selection_and_rank() {
    std::cout << "Running Pivot Selection and Rank tests..." << std::endl;

    linalg::matrix Z(2, 2);
    Z(0,0) = 0.0; Z(0,1) = 1.0;
    Z(1,0) = 0.0; Z(1,1) = 2.0;
    assert(Z.rank() == 1);

    linalg::matrix T(4, 2);
    T(0,0) = 1.0; T(0,1) = 2.0;
    T(1,0) = 2.0; T(1,1) = 4.0;
    T(2,0) = 3.0; T(2,1) = 6.0;
    T(3,0) = 1.0; T(3,1) = -1.0;
    assert(T.rank() == 2);

    linalg::matrix W(2, 4);
    W(0,0) = 1.0; W(0,1) = 2.0; W(0,2) = 3.0; W(0,3) = 4.0;
    W(1,0) = 2.0; W(1,1) = 4.0; W(1,2) = 6.0; W(1,3) = 8.0;
    assert(W.rank() == 1);

    std::cout << "Pivot Selection and Rank tests passed!" << std::endl;
}

void test_precision_row_scaling() {
    std::cout << "Running precision row scaling regression tests..." << std::endl;

    linalg::matrix A(2, 2);
    A(0,0) = 1e-7; A(0,1) = 1.0;
    A(1,0) = 1.0;  A(1,1) = 1.0;

    linalg::matrix invA   = A.inverse();
    linalg::matrix Icheck = A * invA;

    // The condition number of this matrix is ~2e7, so A*A^-1 error is
    // O(machine_eps * cond) ≈ 2e-9; use 1e-6 as a safe, realistic bound.
    assert(is_close(Icheck(0,0), 1.0, 1e-6));
    assert(is_close(Icheck(1,1), 1.0, 1e-6));
    assert(is_close(Icheck(0,1), 0.0, 1e-6));
    assert(is_close(Icheck(1,0), 0.0, 1e-6));

    std::cout << "Precision row scaling tests passed!" << std::endl;
}

void test_rectangular_geometric_utilities() {
    std::cout << "Running rectangular geometry inner product and angle tests..." << std::endl;

    linalg::matrix A(2, 4);
    A(0,0) = 1.0;

    linalg::matrix B(2, 4);
    B(0,1) = 1.0;

    assert(is_close(linalg::inner_product(A, B), 0.0));
    assert(is_close(linalg::angle(A, B), linalg::PI / 2.0));

    linalg::matrix mis(3, 3);
    bool caught = false;
    try { linalg::inner_product(A, mis); } catch (const std::runtime_error&) { caught = true; }
    assert(caught && "Should throw on shape mismatch.");

    std::cout << "Rectangular geometric utility tests passed!" << std::endl;
}

void test_to_string_formatting() {
    std::cout << "Running to_string formatting tests..." << std::endl;

    linalg::matrix A(2, 2);
    A(0,0) = 1.0; A(0,1) = 2.0;
    A(1,0) = 3.0; A(1,1) = 4.0;

    std::string out      = A.to_string();
    std::string expected = "1 2\n3 4\n";
    assert(out == expected && "to_string output did not match expected formatting!");

    std::cout << "to_string formatting tests passed!" << std::endl;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
    std::cout << "=== STARTING ELDA CORE UNIT TESTS ===" << std::endl;

    // Issue #6 – tolerant comparison tests
    test_approx_equal_helper();
    test_matrices_approx_equal();
    test_check_ortho_rotation();
    test_check_ortho_perturbed();
    test_triangular_checks_tolerance();

    // Regression / existing tests
    test_matrix_construction_and_identity();
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

    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}