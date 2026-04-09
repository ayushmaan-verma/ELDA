#ifndef ELDA_TRANSFORMS_H
#define ELDA_TRANSFORMS_H

#include "matrix.hpp"

namespace linalg {

    /// Returns a 3x3 homogeneous translation matrix for 2D coordinates.
    matrix shift(double dx, double dy);
    /// Returns a 3x3 homogeneous scaling matrix for 2D coordinates.
    matrix scale(double kx, double ky);
    /// Returns a 3x3 homogeneous rotation matrix for 2D coordinates.
    /// `angle` is measured in radians.
    matrix rotate(double angle);

    /// Returns a 4x4 homogeneous translation matrix for 3D coordinates.
    matrix shift(double dx, double dy, double dz);
    /// Returns a 4x4 homogeneous scaling matrix for 3D coordinates.
    matrix scale(double kx, double ky, double kz);
    /// Returns a 4x4 homogeneous rotation matrix around the x-axis.
    /// `angle` is measured in radians.
    matrix rot_x(double angle);
    /// Returns a 4x4 homogeneous rotation matrix around the y-axis.
    /// `angle` is measured in radians.
    matrix rot_y(double angle);
    /// Returns a 4x4 homogeneous rotation matrix around the z-axis.
    /// `angle` is measured in radians.
    matrix rot_z(double angle);

}
#endif // ELDA_TRANSFORMS_H
