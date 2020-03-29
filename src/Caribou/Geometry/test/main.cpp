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

template <typename Derived>
::testing::AssertionResult CmpHelperMatrixEQ(const char* lhs_expression,
                                             const char* rhs_expression,
                                             const char* /*abs_error_expr*/,
                                             const Eigen::MatrixBase<Derived> & lhs_value,
                                             const Eigen::MatrixBase<Derived> & rhs_value,
                                             double abs_error) {
    auto max_error = ((lhs_value) - (rhs_value)).cwiseAbs().maxCoeff();

    if (max_error < abs_error) {
        return ::testing::AssertionSuccess();
    }

    Eigen::IOFormat clean(static_cast<int>(std::abs(std::log10(abs_error))), 0, " ", "\n", "    |", "|");

    return ::testing::AssertionFailure()
        << "The difference between " << lhs_expression << " and " << rhs_expression
        << " is " << max_error << ", which exceeds " << abs_error << ", where\n"
        << lhs_expression << " evaluates to \n" << lhs_value.format(clean) << "\n"
        << rhs_expression << " evaluates to \n" << rhs_value.format(clean) << "\n";
}

#define EXPECT_MATRIX_NEAR(val1, val2, abs_error)\
  EXPECT_PRED_FORMAT3(::CmpHelperMatrixEQ, \
                      val1, val2, abs_error)

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
