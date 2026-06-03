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
    assert(is_close(m.arr[0][0], 0.0));

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
    assert(is_close(sum.arr[0][0], 5.0));
    assert(is_close(sum.arr[1][1], 5.0));
    linalg::matrix scaled = m1 * 2.0;
    assert(is_close(scaled.arr[0][0], 4.0));
    std::cout << "Arithmetic tests passed!" << std::endl;
}

void test_matrix_transpose_and_determinant() {
    std::cout << "Running transpose & determinant tests..." << std::endl;
    linalg::matrix m(2, 2);
    m.arr = {{1.0, 2.0}, {3.0, 4.0}};
    assert(is_close(m.det(), -2.0));
    linalg::matrix t = m.transpose();
    assert(is_close(t.arr[0][1], 3.0));
    assert(is_close(t.arr[1][0], 2.0));
    std::cout << "Transpose & determinant tests passed!" << std::endl;
}

void test_transforms_and_vectors() {
    std::cout << "Running transforms & vector utilities tests..." << std::endl;
    linalg::matrix trans = linalg::shift(5.0, 10.0);
    assert(is_close(trans.arr[0][2], 5.0));
    assert(is_close(trans.arr[1][2], 10.0));
    linalg::matrix v3 = linalg::vec3(4.5, 5.5, 6.5);
    assert(v3.row == 3);
    assert(v3.col == 1);
    assert(is_close(v3.arr[0][0], 4.5));
    assert(is_close(v3.arr[2][0], 6.5));
    std::cout << "Transforms & vector utilities tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA UNIT TESTS ===" << std::endl;
    test_matrix_construction_and_identity();
    test_matrix_arithmetic();
    test_matrix_transpose_and_determinant();
    test_transforms_and_vectors();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
