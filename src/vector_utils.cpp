#include "elda/vector_utils.hpp"

namespace linalg {

// Convenience constructors for zero-initialized column vectors.
matrix vec1() { return matrix(1, 1); }

matrix vec2() { return matrix(2, 1); }

matrix vec3() { return matrix(3, 1); }

matrix vec4() { return matrix(4, 1); }

matrix vec5() { return matrix(5, 1); }

matrix vec1(double x) {
    matrix m(1, 1);
    m.arr[0][0] = x;
    return m;
}

matrix vec2(double x, double y) {
    matrix m(2, 1);
    m.arr[0][0] = x;
    m.arr[1][0] = y;
    return m;
}

matrix vec3(double x, double y, double z) {
    matrix m(3, 1);
    m.arr[0][0] = x;
    m.arr[1][0] = y;
    m.arr[2][0] = z;
    return m;
}

matrix vec4(double x, double y, double z, double w) {
    matrix m(4, 1);
    m.arr[0][0] = x;
    m.arr[1][0] = y;
    m.arr[2][0] = z;
    m.arr[3][0] = w;
    return m;
}

matrix vec5(double x, double y, double z, double w, double u) {
    matrix m(5, 1);
    m.arr[0][0] = x;
    m.arr[1][0] = y;
    m.arr[2][0] = z;
    m.arr[3][0] = w;
    m.arr[4][0] = u;
    return m;
}

bool check_lin_comb(matrix m1, matrix m2) {
    // m2 is in span{m1} iff rank([m1 m2]) == rank([m1]).
    matrix aug(m1.row, 2);
    for (int i = 0; i < m1.row; i++) {
        aug.arr[i][0] = m1.arr[i][0];
        aug.arr[i][1] = m2.arr[i][0];
    }

    matrix coeff(m1.row, 1);
    for (int i = 0; i < m1.row; i++) {
        coeff.arr[i][0] = m1.arr[i][0];
    }

    return aug.rank() == coeff.rank();
}

bool check_lin_comb(matrix m1, matrix m2, matrix m3) {
    // m3 is in span{m1, m2} iff rank([m1 m2 m3]) == rank([m1 m2]).
    matrix aug(m1.row, 3);
    for (int i = 0; i < m1.row; i++) {
        aug.arr[i][0] = m1.arr[i][0];
        aug.arr[i][1] = m2.arr[i][0];
        aug.arr[i][2] = m3.arr[i][0];
    }

    matrix coeff(m1.row, 2);
    for (int i = 0; i < m1.row; i++) {
        coeff.arr[i][0] = m1.arr[i][0];
        coeff.arr[i][1] = m2.arr[i][0];
    }

    return aug.rank() == coeff.rank();
}

bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4) {
    // m4 is in span{m1, m2, m3} iff rank([m1 m2 m3 m4]) == rank([m1 m2 m3]).
    matrix aug(m1.row, 4);
    for (int i = 0; i < m1.row; i++) {
        aug.arr[i][0] = m1.arr[i][0];
        aug.arr[i][1] = m2.arr[i][0];
        aug.arr[i][2] = m3.arr[i][0];
        aug.arr[i][3] = m4.arr[i][0];
    }

    matrix coeff(m1.row, 3);
    for (int i = 0; i < m1.row; i++) {
        coeff.arr[i][0] = m1.arr[i][0];
        coeff.arr[i][1] = m2.arr[i][0];
        coeff.arr[i][2] = m3.arr[i][0];
    }

    return aug.rank() == coeff.rank();
}

bool check_lin_comb(matrix m1, matrix m2, matrix m3, matrix m4, matrix m5) {
    // m5 is in span{m1, m2, m3, m4} iff rank([m1 m2 m3 m4 m5]) == rank([m1 m2 m3 m4]).
    matrix aug(m1.row, 5);
    for (int i = 0; i < m1.row; i++) {
        aug.arr[i][0] = m1.arr[i][0];
        aug.arr[i][1] = m2.arr[i][0];
        aug.arr[i][2] = m3.arr[i][0];
        aug.arr[i][3] = m4.arr[i][0];
        aug.arr[i][4] = m5.arr[i][0];
    }

    matrix coeff(m1.row, 4);
    for (int i = 0; i < m1.row; i++) {
        coeff.arr[i][0] = m1.arr[i][0];
        coeff.arr[i][1] = m2.arr[i][0];
        coeff.arr[i][2] = m3.arr[i][0];
        coeff.arr[i][3] = m4.arr[i][0];
    }

    return aug.rank() == coeff.rank();
}
}
