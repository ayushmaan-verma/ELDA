#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <elda/linalg.hpp>

using linalg::matrix;

bool is_close(double a, double b, double epsilon = 1e-5) {
    return std::abs(a - b) < epsilon;
}

matrix make_augmented(std::vector<std::vector<double>> values) {
    matrix m(values.size(), values[0].size());
    m.arr = values;
    return m;
}

void test_matrix_construction_and_identity() {
    std::cout << "Running construction tests..." << std::endl;
    linalg::matrix m(2, 2);
    assert(m.row == 2);
    assert(m.col == 2);
    assert(is_close(m.arr[0][0], 0.0));
    linalg::matrix id = linalg::identity(3);
    assert(is_close(id.arr[0][0], 1.0));
    assert(is_close(id.arr[1][1], 1.0));
    assert(is_close(id.arr[2][2], 1.0));
    std::cout << "Construction tests passed!" << std::endl;
}

void test_matrix_arithmetic() {
    std::cout << "Running arithmetic tests..." << std::endl;
    linalg::matrix m1(2, 2), m2(2, 2);
    m1.arr = {{1.0, 2.0}, {3.0, 4.0}};
    m2.arr = {{1.0, 1.0}, {1.0, 1.0}};
    linalg::matrix sum = m1 + m2;
    assert(is_close(sum.arr[0][0], 2.0));
    assert(is_close(sum.arr[1][1], 5.0));
    linalg::matrix scaled = m1 * 2.0;
    assert(is_close(scaled.arr[0][0], 2.0));
    assert(is_close(scaled.arr[1][1], 8.0));
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

void test_transforms() {
    std::cout << "Running transforms tests..." << std::endl;
    linalg::matrix trans = linalg::shift(5.0, 10.0);
    assert(trans.row == 3);
    assert(trans.col == 3);
    assert(is_close(trans.arr[0][2], 5.0));
    assert(is_close(trans.arr[1][2], 10.0));
    linalg::matrix rot = linalg::rotate(0.0);
    assert(is_close(rot.arr[0][0], 1.0));
    assert(is_close(rot.arr[1][1], 1.0));
    std::cout << "Transforms tests passed!" << std::endl;
}

void test_solve_unique_system() {
    std::cout << "Running solve tests (unique system)..." << std::endl;
    // x + 2y = 3 ; 3x + 4y = 5  =>  x = -1, y = 2
    matrix aug = make_augmented({{1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}});
    matrix sol = aug.solve();
    assert(is_close(sol.arr[0][0], -1.0));
    assert(is_close(sol.arr[1][0], 2.0));
    std::cout << "Solve (unique system) passed!" << std::endl;
}

void test_solve_row_swap() {
    std::cout << "Running solve tests (row swap)..." << std::endl;
    // 0x + 1y = 2 ; 1x + 0y = 3  =>  x = 3, y = 2
    matrix aug = make_augmented({{0.0, 1.0, 2.0}, {1.0, 0.0, 3.0}});
    matrix sol = aug.solve();
    assert(is_close(sol.arr[0][0], 3.0));
    assert(is_close(sol.arr[1][0], 2.0));
    std::cout << "Solve (row swap) passed!" << std::endl;
}

void test_solve_inconsistent_system() {
    std::cout << "Running solve tests (inconsistent system)..." << std::endl;
    // x + y = 1 ; 2x + 2y = 5  =>  no solution
    matrix aug = make_augmented({{1.0, 1.0, 1.0}, {2.0, 2.0, 5.0}});
    bool threw = false;
    try {
        aug.solve();
    } catch (const std::runtime_error& e) {
        threw = true;
        assert(std::string(e.what()).find("inconsistent") != std::string::npos);
    }
    assert(threw);
    std::cout << "Solve (inconsistent system) passed!" << std::endl;
}

void test_solve_underdetermined_system() {
    std::cout << "Running solve tests (underdetermined system)..." << std::endl;
    // x + y = 1 ; 2x + 2y = 2  =>  infinitely many solutions
    matrix aug = make_augmented({{1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}});
    bool threw = false;
    try {
        aug.solve();
    } catch (const std::runtime_error& e) {
        threw = true;
        assert(std::string(e.what()).find("singular") != std::string::npos);
    }
    assert(threw);
    std::cout << "Solve (underdetermined system) passed!" << std::endl;
}

void test_solve_shape_check() {
    std::cout << "Running solve tests (shape check)..." << std::endl;
    matrix square(2, 2);
    bool threw = false;
    try {
        square.solve();
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw);
    std::cout << "Solve (shape check) passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING ELDA UNIT TESTS ===" << std::endl;
    test_matrix_construction_and_identity();
    test_matrix_arithmetic();
    test_matrix_transpose_and_determinant();
    test_transforms();
    test_solve_unique_system();
    test_solve_row_swap();
    test_solve_inconsistent_system();
    test_solve_underdetermined_system();
    test_solve_shape_check();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}
