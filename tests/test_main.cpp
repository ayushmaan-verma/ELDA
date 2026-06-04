#include <iostream>
#include <cassert>
#include <cmath>
#include <tuple>
#include "elda/matrix.hpp"

bool is_close(double a, double b, double epsilon = 1e-5) {
    return std::abs(a - b) < epsilon;
}

void test_lu_decomposition() {
    std::cout << "Running LU decomposition tests..." << std::endl;
    
    // Matrix configuration that strictly requires pivot row permutations
    linalg::matrix A(3, 3);
    A.arr[0][0] = 0.0; A.arr[0][1] = 2.0; A.arr[0][2] = 1.0;
    A.arr[1][0] = 1.0; A.arr[1][1] = -2.0; A.arr[1][2] = -3.0;
    A.arr[2][0] = -1.0; A.arr[2][1] = 1.0; A.arr[2][2] = 2.0;

    auto [P, L, U] = A.lu_decomposition();
    
    linalg::matrix LU = L * U;
    linalg::matrix PA = P * A;
    
    // Verify permutation relation: P * A == L * U
    for(size_t i = 0; i < 3; i++) {
        for(size_t j = 0; j < 3; j++) {
            assert(is_close(PA.arr[i][j], LU.arr[i][j]));
        }
    }

    // Verify boundary case handling for rectangular matrices
    linalg::matrix Rect(2, 3);
    bool caught_rect = false;
    try {
        Rect.lu_decomposition();
    } catch (const std::runtime_error&) {
        caught_rect = true;
    }
    assert(caught_rect && "Should throw on rectangular configurations.");

    std::cout << "LU decomposition tests passed!" << std::endl;
}

int main() {
    std::cout << "=== STARTING LU DECOMPOSITION UNIT TESTS ===" << std::endl;
    test_lu_decomposition();
    std::cout << "=== ALL TESTS PASSED SUCCESSFULLY ===" << std::endl;
    return 0;
}