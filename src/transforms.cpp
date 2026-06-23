#include "elda/transforms.hpp"

namespace linalg {

// cppcheck-suppress unusedFunction
matrix shift(double dx, double dy) {
    matrix m = identity(3);
    *m.ref_element(0, 2) = dx;
    *m.ref_element(1, 2) = dy;
    return m;
}

// cppcheck-suppress unusedFunction
matrix scale(double kx, double ky) {
    matrix m = identity(3);
    *m.ref_element(0, 0) = kx;
    *m.ref_element(1, 1) = ky;
    return m;
}

// cppcheck-suppress unusedFunction
matrix rotate(double angle) {
    matrix m = identity(3);
    *m.ref_element(0, 0) = cos(angle);
    *m.ref_element(0, 1) = -sin(angle);
    *m.ref_element(1, 0) = sin(angle);
    *m.ref_element(1, 1) = cos(angle);
    return m;
}

// cppcheck-suppress unusedFunction
matrix shift(double dx, double dy, double dz) {
    matrix m = identity(4);
    *m.ref_element(0, 3) = dx;
    *m.ref_element(1, 3) = dy;
    *m.ref_element(2, 3) = dz;
    return m;
}

// cppcheck-suppress unusedFunction
matrix scale(double kx, double ky, double kz) {
    matrix m = identity(4);
    *m.ref_element(0, 0) = kx;
    *m.ref_element(1, 1) = ky;
    *m.ref_element(2, 2) = kz;
    return m;
}

// cppcheck-suppress unusedFunction
matrix rot_x(double angle) {
    matrix m = identity(4);
    *m.ref_element(1, 1) = cos(angle);
    *m.ref_element(1, 2) = -sin(angle);
    *m.ref_element(2, 1) = sin(angle);
    *m.ref_element(2, 2) = cos(angle);
    return m;
}

// cppcheck-suppress unusedFunction
matrix rot_y(double angle) {
    matrix m = identity(4);
    *m.ref_element(0, 0) = cos(angle);
    *m.ref_element(0, 2) = -sin(angle);
    *m.ref_element(2, 0) = sin(angle);
    *m.ref_element(2, 2) = cos(angle);
    return m;
}

// cppcheck-suppress unusedFunction
matrix rot_z(double angle) {
    matrix m = identity(4);
    *m.ref_element(0, 0) = cos(angle);
    *m.ref_element(0, 1) = -sin(angle);
    *m.ref_element(1, 0) = sin(angle);
    *m.ref_element(1, 1) = cos(angle);
    return m;
}

}  // namespace linalg
