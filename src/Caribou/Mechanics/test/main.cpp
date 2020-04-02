#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

template<int nRows, int nColumns, int Options=0>
using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;


TEST(Mechanics, Strain) {
    using namespace caribou::mechanics;

    using Hexahedron = caribou::geometry::Hexahedron<caribou::Linear>;
    using Mat33 = Matrix<3,3>;
    const auto I = Mat33::Identity();

    Hexahedron h;
    for (const auto & gauss_node : h.gauss_nodes()) {
        const auto x = gauss_node.position;
        const auto dN_dx = h.dL(x);
        const auto U = Matrix<caribou::geometry::traits<Hexahedron>::NumberOfNodesAtCompileTime, 3>::Zero().eval();
        const Mat33 F = elasticity::strain::F(dN_dx, U);
        ASSERT_EQ(I, F);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
