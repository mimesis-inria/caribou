#include <gtest/gtest.h>
#include <iostream>
#include <Caribou/Geometry/Point.h>

TEST(Entity, data) {

}

TEST(Point, constructors) {
    caribou::geometry::Point3D<> p1(1, 2, 3);
    caribou::geometry::Point3D<> p2 = {1, 2, 3};
    caribou::geometry::Point3D<> p3({1, 2, 3});
    caribou::geometry::Point3D<> p4; p4[0] = p1[0]; p4[1]= p2[1]; p4[2] = p3[2];
    caribou::geometry::Point3D<> p5 = p4;

    ASSERT_EQ(p5[0], 1);
    ASSERT_EQ(p5[1], 2);
    ASSERT_EQ(p5[2], 3);
    ASSERT_EQ(p1, p2);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
