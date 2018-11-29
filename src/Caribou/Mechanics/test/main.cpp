#include <gtest/gtest.h>

#include <Caribou/Algebra/Matrix.h>
#include <Caribou/Geometry/LinearHexahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

TEST(Mechanics, Strain) {
    using namespace caribou::geometry;
    using namespace caribou::algebra;
    using namespace caribou::mechanics;

    using Mat33 = Matrix<3,3,FLOATING_POINT_TYPE>;

    const LinearHexahedron initial_hexahedron;
    {
        const Mat33 F = elasticity::strain::F(initial_hexahedron, initial_hexahedron, 0, 0, 0);
        ASSERT_EQ(Mat33::Identity(), F);
    }

    LinearHexahedron deformed_hexahedron = initial_hexahedron;
    deformed_hexahedron.node(8).translate({1,1,1});
    {
        const Mat33 F = elasticity::strain::F(initial_hexahedron, deformed_hexahedron, 0, 0, 0);
        ASSERT_EQ(Mat33 ({{0.7999999523, -0.2000000179, -0.2000000179}, {-0.2000000179, 0.7999999523, -0.2000000179}, {-0.200000003, -0.200000003, 0.8000000119}}), F);
    }

}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
