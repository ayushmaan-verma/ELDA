#include <elda/linalg.hpp>
#include <iostream>

using namespace linalg;

int main() {
    matrix m(2, 2, 5.0);
    std::cout << m;
}