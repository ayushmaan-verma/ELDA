#ifndef ELDA_VECTOR_UTILS_H
#define ELDA_VECTOR_UTILS_H

#include "matrix.hpp"

namespace linalg {

    /// Returns a 1D column vector.
    matrix vec1(double x);
    /// Returns a 2D column vector.
    matrix vec2(double x, double y);
    /// Returns a 3D column vector.
    matrix vec3(double x, double y, double z);
    /// Returns a 4D column vector.
    matrix vec4(double x, double y, double z, double w);
    /// Returns a 5D column vector.
    matrix vec5(double x, double y, double z, double w, double u);

}

#endif // ELDA_VECTOR_UTILS_H
