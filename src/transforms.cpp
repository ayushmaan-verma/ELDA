#include "elda/transforms.hpp"

namespace linalg {

    matrix shift(double dx, double dy) {
        matrix m = identity(3);
        *m.ref_element(0,2) = dx;
        *m.ref_element(1,2) = dy;
        return m;
    }
    matrix scale(double kx, double ky) {
        matrix m = identity(3);
        *m.ref_element(0,0) = kx;
        *m.ref_element(1,1) = ky;
        return m;
    }
    matrix rotate(double angle) {
        matrix m = identity(3);
        *m.ref_element(0,0) = cos(angle);
        *m.ref_element(0,1) = -sin(angle);
        *m.ref_element(1,0) = sin(angle);
        *m.ref_element(1,1) = cos(angle);
        return m;
    }

    matrix shift(double dx, double dy, double dz) {
        matrix m = identity(4);
        *m.ref_element(0,3) = dx;
        *m.ref_element(1,3) = dy;
        *m.ref_element(2,3) = dz;
        return m;
    }
    matrix scale(double kx, double ky, double kz) {
        matrix m = identity(4);
        *m.ref_element(0,0) = kx;
        *m.ref_element(1,1) = ky;
        *m.ref_element(2,2) = kz;
        return m;
    }
    matrix rot_x(double angle) {
        matrix m = identity(4);
        *m.ref_element(1,1) = cos(angle);
        *m.ref_element(1,2) = -sin(angle);
        *m.ref_element(2,1) = sin(angle);
        *m.ref_element(2,2) = cos(angle);
        return m;
    }
    matrix rot_y(double angle) {
        matrix m = identity(4);
        *m.ref_element(0,0) = cos(angle);
        *m.ref_element(0,2) = -sin(angle);
        *m.ref_element(2,0) = sin(angle);
        *m.ref_element(2,2) = cos(angle);
        return m;
    }
    matrix rot_z(double angle) {
        matrix m = identity(4);
        *m.ref_element(0,0) = cos(angle);
        *m.ref_element(0,1) = -sin(angle);
        *m.ref_element(1,0) = sin(angle);
        *m.ref_element(1,1) = cos(angle);
        return m;
    }

}
