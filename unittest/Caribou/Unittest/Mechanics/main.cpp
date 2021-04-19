#include <gtest/gtest.h>
#include <Caribou/config.h>
#include <Eigen/Core>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

template <typename Derived>
::testing::AssertionResult CmpHelperMatrixEQ(const char* lhs_expression,
    const char* rhs_expression,
    const char* /*abs_error_expr*/,
    const Eigen::MatrixBase<Derived>& lhs_value,
    const Eigen::MatrixBase<Derived>& rhs_value,
    double abs_error) {
    auto max_error = ((lhs_value)-(rhs_value)).cwiseAbs().maxCoeff();

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


TEST(Mechanics, Strain) {
    using namespace caribou::mechanics;

    using Hexahedron = caribou::geometry::Hexahedron<caribou::Linear>;
    using Mat33 = Eigen::Matrix<FLOATING_POINT_TYPE, 3,3>;
    const Mat33 I = Mat33::Identity();

    Hexahedron h;
    for (const auto & gauss_node : h.gauss_nodes()) {
        const auto x = gauss_node.position;
        const auto dN_dx = h.dL(x);
        const auto U = Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3>::Zero().eval();
        const Mat33 F = elasticity::strain::F(dN_dx, U);
        EXPECT_MATRIX_NEAR(I, F, 1e-15);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
