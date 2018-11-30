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
    deformed_hexahedron.node(7).translate({1,1,1});
    {
        const Mat33 F = elasticity::strain::F(initial_hexahedron, deformed_hexahedron, 0, 0, 0);
        ASSERT_EQ(Mat33 ({{1.125, 0.125, 0.125}, {-0.125, 0.875, -0.125}, {-0.125, -0.125, 0.875}}), F);

        const Mat33 e = elasticity::strain::small_strain(initial_hexahedron, deformed_hexahedron, LinearHexahedron().node(7));
        std::cout << e;
        std::cout << Vector<3>({1,1,1}).unit() << "\n";
//        std::cout << (e * Vector<3>({1,1,1}).unit()) * Vector<3>({1,1,1}).unit();
    }

}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
