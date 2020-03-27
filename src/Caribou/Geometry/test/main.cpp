#include <gtest/gtest.h>
#include <Caribou/config.h>
#include <Eigen/Core>

template <class Derived>
struct is_eigen : public std::is_base_of<Eigen::DenseBase<Derived>, Derived> {
};
template <class Derived,
    class = typename std::enable_if<is_eigen<Derived>::value>::type>
::std::ostream &operator<<(::std::ostream &o, const Derived &m) {
    o << "\n" << static_cast<const Eigen::DenseBase<Derived> &>(m);
    return o;
}

template <typename LocalCoordinates>
FLOATING_POINT_TYPE p1(const LocalCoordinates & p) {
    if (p.rows() == 1)  return 5 + 2*p[0];
    if (p.rows() == 2)  return 5 + 2*p[0] + 3*p[1];
    if (p.rows() == 3)  return 5 + 2*p[0] + 3*p[1] + 4*p[2];
    return 0;
}

template <typename LocalCoordinates>
FLOATING_POINT_TYPE p2(const LocalCoordinates & p) {
    if (p.rows() == 1)  return 5 + 2*p[0]*p[0];
    if (p.rows() == 2)  return 5 + 2*p[0]*p[1] + 3*p[1]*p[1];
    if (p.rows() == 3)  return 5 + 2*p[0]*p[1] + 3*p[1]*p[2]+ 4*p[2]*p[2];
    return 0;
}

#include "test_segment.h"
#include "test_triangle.h"
#include "test_quad.h"
#include "test_tetrahedron.h"
#include "test_hexahedron.h"

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
