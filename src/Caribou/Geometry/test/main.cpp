#include <gtest/gtest.h>
#include <Caribou/Geometry/Point.h>
#include <Caribou/Geometry/Segment.h>

struct Data {
    char a;
    char b;
    Data(): a(0), b(0) {}
    Data(char a, char b): a(a), b(b) {}
    bool operator==(const Data & d) const {
        return (
                this->a == d.a && this->b == d.b
        );
    }
};

TEST(Geometry, Point) {
    auto p1 = caribou::geometry::make_point({1, 2, 3});
    auto p2 = caribou::geometry::make_point(1, 2, 3);

    ASSERT_EQ(p1, p2);
    ASSERT_EQ(p1*p2, 14);

    caribou::geometry::Point3D<caribou::geometry::BaseData, int> p3({1, 2, 3});
    caribou::geometry::Point3D<caribou::geometry::BaseData, int> p4;
    p4[0] = p1[0]; p4[1]= p2[1]; p4[2] = p3[2];

    ASSERT_EQ(p3, p4);
    ASSERT_EQ(p1, p4);

    auto p5 = p4;

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

    // Make sure a point have minimal space taken in the memory
    ASSERT_EQ(sizeof(p8), (size_t) 5);
}

TEST(Geometry, Segment) {
    using namespace caribou::geometry;
    using Point3D = Point3D<>;

    Point3D p1(0, 0, 0);
    Point3D p2 = {1, 1, 1};
    auto s1 =  make_segment(p1, p2);

    Point3D p3(0, 0, 0);
    Point3D p4 = {1, 1, 1};
    auto s2 = make_segment(&p3, &p4);

    ASSERT_EQ(s1, s2);
    ASSERT_EQ(s2, s1);
    p3 = {9, 9, 9};

    ASSERT_NE(s1, s2);
    ASSERT_NE(s2, s1);

    // Make sure a segment have minimal space taken in the memory
    caribou::geometry::Point3D<Data, char> p5;
    caribou::geometry::Point3D<Data, char> p6;
    Data d = {(char) 1, (char)2};
    auto s3 = make_segment(p5, p6, d);
    ASSERT_EQ(
            sizeof(s3),

            3*sizeof(char) /* p5 3xcoordinates */ + 2*sizeof(char) /* p5 data */ +
            3*sizeof(char) /* p6 3xcoordinates */ + 2*sizeof(char) /* p6 data */ +
            2*sizeof(char) /* segment data */
    );
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
