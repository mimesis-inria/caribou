#pragma once

#include <Eigen/Core>

namespace caribou::algebra {

/**
 * Given a second order symmetric tensor a (3x3 matrix), gives the result of
 *      M_ijkl = a_ij a_kl
 * which is a super-symmetric 4th order tensor (major and minor symmetries) represented by a 6x6 matrix.
 *
 * \warning Only works if a is symmetric (a_ij == a_ji)
 */
template <typename Real, int options>
auto symmetric_dyad_1 (const Eigen::Matrix<Real, 3, 3, options> & a) {
    Eigen::Matrix<Real, 6, 6, options> M;
    M <<
        a(0,0)*a(0,0),    a(0,0)*a(1,1),    a(0,0)*a(2,2),    a(0,0)*a(0,1),    a(0,0)*a(1,2),    a(0,0)*a(0,2),
        a(0,0)*a(1,1),    a(1,1)*a(1,1),    a(1,1)*a(2,2),    a(1,1)*a(0,1),    a(1,1)*a(1,2),    a(1,1)*a(0,2),
        a(0,0)*a(2,2),    a(1,1)*a(2,2),    a(2,2)*a(2,2),    a(2,2)*a(0,1),    a(2,2)*a(1,2),    a(2,2)*a(0,2),
        a(0,0)*a(0,1),    a(1,1)*a(0,1),    a(2,2)*a(0,1),    a(0,1)*a(0,1),    a(0,1)*a(1,2),    a(0,1)*a(0,2),
        a(0,0)*a(1,2),    a(1,1)*a(1,2),    a(2,2)*a(1,2),    a(0,1)*a(1,2),    a(1,2)*a(1,2),    a(1,2)*a(0,2),
        a(0,0)*a(0,2),    a(1,1)*a(0,2),    a(2,2)*a(0,2),    a(0,1)*a(0,2),    a(1,2)*a(0,2),    a(0,2)*a(0,2);
    return M;
}

/**
 * Given a second order symmetric tensor a (3x3 matrix), gives the result of
 *      M_ijkl = 1/2 * (a_ik a_jl + a_il a_jk)
 * which is a super-symmetric 4th order tensor (major and minor symmetries) represented by a 6x6 matrix.
 *
 * \warning Only works if a is symmetric (a_ij == a_ji)
 */
template <typename Real, int options>
auto symmetric_dyad_2 (const Eigen::Matrix<Real, 3, 3, options> & a) {
    Eigen::Matrix<Real, 6, 6, options> M;
    M <<
        a(0,0)*a(0,0),    a(0,1)*a(0,1),    a(0,2)*a(0,2),                    a(0,0)*a(0,1),                             a(0,1)*a(0,2),                             a(0,0)*a(0,2),
        a(0,1)*a(0,1),    a(1,1)*a(1,1),    a(1,2)*a(1,2),                    a(1,1)*a(0,1),                             a(1,1)*a(1,2),                             a(1,2)*a(0,1),
        a(0,2)*a(0,2),    a(1,2)*a(1,2),    a(2,2)*a(2,2),                    a(1,2)*a(0,2),                             a(2,2)*a(1,2),                             a(2,2)*a(0,2),
        a(0,0)*a(0,1),    a(1,1)*a(0,1),    a(1,2)*a(0,2),    1/2. * (a(0,0)*a(1,1) + a(0,1)*a(0,1)),    1/2. * (a(0,1)*a(1,2) + a(0,2)*a(1,1)),    1/2. * (a(0,0)*a(1,2) + a(0,2)*a(0,1)),
        a(0,1)*a(0,2),    a(1,1)*a(1,2),    a(2,2)*a(1,2),    1/2. * (a(0,1)*a(1,2) + a(0,2)*a(1,1)),    1/2. * (a(1,1)*a(2,2) + a(1,2)*a(1,2)),    1/2. * (a(0,1)*a(2,2) + a(1,2)*a(0,2)),
        a(0,0)*a(0,2),    a(1,2)*a(0,1),    a(2,2)*a(0,2),    1/2. * (a(0,0)*a(1,2) + a(0,2)*a(0,1)),    1/2. * (a(0,1)*a(2,2) + a(1,2)*a(0,2)),    1/2. * (a(0,0)*a(2,2) + a(0,2)*a(0,2));
    return M;
}

} // namespace caribou::algebra