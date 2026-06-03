#include <iostream>
#include <cassert>
#include <cmath>
#include <elda/linalg.hpp>

bool is_close(double a, double b, double epsilon = 1e-5) {
    return std::abs(a - b) < epsilon;
}

void test_matrix_construction_and_identity() {
    std::cout << "Running construction tests..." << std::endl;
    linalg::matrix m(2, 2, 0.0);
    assert(m.get_rows() == 2);
    assert(m.get_cols() == 2);
    assert(is_close(m(0, 0), 0.0));
    std::cout << "Construction tests passed!" << std::endl;
}

void test_matrix_arithmetic() {
    std::cout << "Running arithmetic tests..." << std::endl;
    linalg::matrix m1(2, 2, 2.0);
    linalg::matrix m2(2, 2, 3.0);
    linalg::matrix sum = m1 + m2;
    assert(is_close(sum(0, 0), 5.0));
    assert(is_close(sum(1, 1), 5.0));
    linalg::matrix scaled = m1 * 2.0;
    assert(is_close(scaled(0, 0), 4.0));
    std::cout << "Arithmetic tests passed!" << std::endl;
}

void test_matrix_transpose_and_determinant() {
    std::cout << "Running transpose & determinant tests..." << std::endl;
    linalg::matrix m(2, 2, 0.0);
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;
    assert(is_close(m.det(), -2.0));
    linalg::matrix t = m.transpose();
    assert(is_close(t(0, 1), 3.0));
    assert(is_close(t(1, 0), 2.0));
    std::cout << "Transpose & determinant tests passed!" << std::endl;
}

void test_transforms_and_vectors() {
    std::cout << "Running transforms & vector utilities tests..." << std::endl;
    linalg::matrix trans = linalg::translation2d(5.0, 10.0);
    assert(is_close(trans(0, 2), 5.0));
    assert(is_close(trans(1, 2), 10.0));
    linalg::matrix v3 = linalg::vec3(4.5);
    assert(v3.get_rows() == 3);
    assert(v3.get_cols() == 1);
    assert(is_close(v3(0, 0), 4.5));
    assert(is_close(v3(2, 0), 4.5));
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

int main() {
    std::cout << "=== STARTING ELDA UNIT TESTS ===" << std::endl;
    test_matrix_construction_and_identity();
    test_matrix_arithmetic();
    test_matrix_transpose_and_determinant();
    test_transforms_and_vectors();
    test_vector_shape_validation();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
