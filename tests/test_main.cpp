#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <elda/linalg.hpp>

bool is_close(double a, double b, double epsilon = 1e-5) {
    return std::abs(a - b) < epsilon;
}

void test_matrix_construction_and_identity() {
    std::cout << "Running construction tests..." << std::endl;
    linalg::matrix m(2, 2);
    assert(m.row == 2);
    assert(m.col == 2);
    assert(is_close(m.get_element(0, 0), 0.0));
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
    linalg::matrix v3 = linalg::vec3(4.5, 4.5, 4.5);
    assert(v3.row == 3);
    assert(v3.col == 1);
    assert(is_close(v3.get_element(0, 0), 4.5));
    assert(is_close(v3.get_element(2, 0), 4.5));
    std::cout << "Transforms & vector utilities tests passed!" << std::endl;
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
    test_matpow_validation_and_results();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
