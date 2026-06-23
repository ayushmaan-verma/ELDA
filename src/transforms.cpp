#include "elda/transforms.hpp"

namespace linalg {

[[maybe_unused]] matrix shift(double dx, double dy) {
    matrix m = identity(3);
    *m.ref_element(0, 2) = dx;
    *m.ref_element(1, 2) = dy;
    return m;
}

[[maybe_unused]] matrix scale(double kx, double ky) {
    matrix m = identity(3);
    *m.ref_element(0, 0) = kx;
    *m.ref_element(1, 1) = ky;
    return m;
}

[[maybe_unused]] matrix rotate(double theta) {
    matrix m = identity(3);
    *m.ref_element(0, 0) = cos(theta);
    *m.ref_element(0, 1) = -sin(theta);
    *m.ref_element(1, 0) = sin(theta);
    *m.ref_element(1, 1) = cos(theta);
    return m;
}

[[maybe_unused]] matrix shift(double dx, double dy, double dz) {
    matrix m = identity(4);
    *m.ref_element(0, 3) = dx;
    *m.ref_element(1, 3) = dy;
    *m.ref_element(2, 3) = dz;
    return m;
}

[[maybe_unused]] matrix scale(double kx, double ky, double kz) {
    matrix m = identity(4);
    *m.ref_element(0, 0) = kx;
    *m.ref_element(1, 1) = ky;
    *m.ref_element(2, 2) = kz;
    return m;
}

[[maybe_unused]] matrix rot_x(double theta) {
    matrix m = identity(4);
    *m.ref_element(1, 1) = cos(theta);
    *m.ref_element(1, 2) = -sin(theta);
    *m.ref_element(2, 1) = sin(theta);
    *m.ref_element(2, 2) = cos(theta);
    return m;
}

[[maybe_unused]] matrix rot_y(double theta) {
    matrix m = identity(4);
    *m.ref_element(0, 0) = cos(theta);
    *m.ref_element(0, 2) = -sin(theta);
    *m.ref_element(2, 0) = sin(theta);
    *m.ref_element(2, 2) = cos(theta);
    return m;
}

[[maybe_unused]] matrix rot_z(double theta) {
    matrix m = identity(4);
    *m.ref_element(0, 0) = cos(theta);
    *m.ref_element(0, 1) = -sin(theta);
    *m.ref_element(1, 0) = sin(theta);
    *m.ref_element(1, 1) = cos(theta);
    return m;
}

}  // namespace linalg
