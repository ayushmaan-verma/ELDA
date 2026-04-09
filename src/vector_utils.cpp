#include "elda/vector_utils.hpp"

namespace linalg {

    matrix vec1(double x) {
        matrix m (1,1);
        m.arr[0][0] = x;
        return m;
    }
    matrix vec2(double x, double y) {
        matrix m (2,1);
        m.arr[0][0] = x;
        m.arr[1][0] = y;
        return m;
    }
    matrix vec3(double x, double y, double z) {
        matrix m (3,1);
        m.arr[0][0] = x;
        m.arr[1][0] = y;
        m.arr[2][0] = z;
        return m;
    }
    matrix vec4(double x, double y, double z, double w) {
        matrix m (4,1);
        m.arr[0][0] = x;
        m.arr[1][0] = y;
        m.arr[2][0] = z;
        m.arr[3][0] = w;
        return m;
    }
    matrix vec5(double x, double y, double z, double w, double u) {
        matrix m (5,1);
        m.arr[0][0] = x;
        m.arr[1][0] = y;
        m.arr[2][0] = z;
        m.arr[3][0] = w;
        m.arr[4][0] = u;
        return m;
    }
}
