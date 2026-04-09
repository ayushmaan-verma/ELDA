#ifndef ELDA_VECTOR_UTILS_H
#define ELDA_VECTOR_UTILS_H

#include "matrix.hpp"

namespace linalg {

    /// Returns a zero-initialized 1D column vector (`R^1`).
    matrix vec1();
    /// Returns a zero-initialized 2D column vector (`R^2`).
    matrix vec2();
    /// Returns a zero-initialized 3D column vector (`R^3`).
    matrix vec3();
    /// Returns a zero-initialized 4D column vector (`R^4`).
    matrix vec4();
    /// Returns a zero-initialized 5D column vector (`R^5`).
    matrix vec5();

    /// Returns a 1D column vector with entry `x`.
    matrix vec1(double x);
    /// Returns a 2D column vector with entries `(x, y)`.
    matrix vec2(double x, double y);
    /// Returns a 3D column vector with entries `(x, y, z)`.
    matrix vec3(double x, double y, double z);
    /// Returns a 4D column vector with entries `(x, y, z, w)`.
    matrix vec4(double x, double y, double z, double w);
    /// Returns a 5D column vector with entries `(x, y, z, w, u)`.
    matrix vec5(double x, double y, double z, double w, double u);

    /// Returns true if `m2` is a linear combination of `m1`.
    bool check_lin_comb(matrix m1, matrix m2);
    /// Returns true if `m3` is a linear combination of `m1` and `m2`.
    bool check_lin_comb(matrix m1, matrix m2, matrix m3);
    /// Returns true if `m4` is a linear combination of `m1`, `m2`, and `m3`.
    bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4);
    /// Returns true if `m5` is a linear combination of `m1`, `m2`, `m3`, and `m4`.
    bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4, matrix m5);
}

#endif // ELDA_VECTOR_UTILS_H
