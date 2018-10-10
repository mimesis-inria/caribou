#include <gtest/gtest.h>
#include <Caribou/Geometry/Point.h>
#include <Caribou/Geometry/Segment.h>

struct Data {
    char a;
    char b;
    bool operator==(const Data & d) const {
        return (
                this->a == d.a && this->b == d.b
        );
    }
};

TEST(Geometry, Point) {
    caribou::geometry::Point3D<> p1(1, 2, 3);
    caribou::geometry::Point3D<> p2 = {1, 2, 3};
    caribou::geometry::Point3D<> p3({1, 2, 3});
    caribou::geometry::Point3D<> p4; p4[0] = p1[0]; p4[1]= p2[1]; p4[2] = p3[2];
    caribou::geometry::Point3D<> p5 = p4;

    ASSERT_EQ(p5.x(), 1);

    p5.set_y(4);
    ASSERT_EQ(p5[1], 4);

    ASSERT_EQ(p5[2], 3);
    ASSERT_EQ(p1, p2);

    caribou::geometry::Point3D<Data, char> p6 (0, 0, 0);
    p6.data.a = (char)1;
    p6.data.b = (char)2;

    caribou::geometry::Point3D<Data, char> p7 (0, 0, 0);
    p7.data.a = (char)1;
    p7.data.b = (char)2;

    caribou::geometry::Point3D<Data, char> p8;
    p8.data.a = 5;
    p8.data.b = 6;

    ASSERT_EQ(p6, p7);
    ASSERT_NE(p6, p8);

    ASSERT_EQ(sizeof(p8), (size_t) 5);
}

TEST(Geometry, Segment) {
    caribou::geometry::Point3D<> p1(0, 0, 0);
    caribou::geometry::Point3D<> p2 = {1, 1, 1};
    caribou::geometry::Segment<> s1 =  { p1, p2 };

    caribou::geometry::Point3D<> p3(0, 0, 0);
    caribou::geometry::Point3D<> p4 = {1, 1, 1};
    caribou::geometry::Segment<> s2 ( p3, p4 );

    ASSERT_EQ(s1, s2);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
