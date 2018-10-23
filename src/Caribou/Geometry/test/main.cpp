#include <gtest/gtest.h>
#include <Caribou/Geometry/Polygon.h>

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

    caribou::geometry::Point3D<caribou::geometry::BaseData, unsigned char> p3({1, 2, 3});
    caribou::geometry::Point3D<caribou::geometry::BaseData, unsigned char> p4;
    p4[0] = p1[0]; p4[1]= p2[1]; p4[2] = p3[2];

    ASSERT_EQ(p3, p4);
    ASSERT_EQ(p1, p4);

    auto p5 = p4;

    ASSERT_EQ(p5.x(), 1);

    p5.set_y(4);
    ASSERT_EQ(p5[1], 4);

    ASSERT_EQ(p5[2], 3);
    ASSERT_EQ(p1, p2);

    caribou::geometry::Point3D<Data, unsigned char> p6 (0, 0, 0);
    p6.data.a = (char)1;
    p6.data.b = (char)2;

    caribou::geometry::Point3D<Data, unsigned char> p7 (0, 0, 0);
    p7.data.a = (char)1;
    p7.data.b = (char)2;

    caribou::geometry::Point3D<Data, unsigned char> p8;
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
    caribou::geometry::Point3D<Data, unsigned char> p5;
    caribou::geometry::Point3D<Data, unsigned char> p6;
    Data d = {(unsigned char) 1, (unsigned char)2};
    auto s3 = make_segment(p5, p6, d);
    ASSERT_EQ(
            sizeof(s3),

            3*sizeof(unsigned char) /* p5 3xcoordinates */ + 2*sizeof(unsigned char) /* p5 data */ +
            3*sizeof(unsigned char) /* p6 3xcoordinates */ + 2*sizeof(unsigned char) /* p6 data */ +
            2*sizeof(unsigned char) /* segment data */
    );
}

TEST(Geometry, Polygon) {
    using namespace caribou::geometry;
    using Point3D = Point3D<>;

    Point3D p1 = {0, 0, 0};
    Point3D p2 = {1, 0, 0};
    Point3D p3 = {0, 1, 0};

    auto s1 = make_segment(p1, p2);
    auto s2 = make_segment(p2, p3);
    auto s3 = make_segment(p3, p1);

    // Create a polygon from 3 segments
    auto polygon_1 = make_polygon(s1, s2, s3);

    // Create a copied polygone
    auto polygon_2 = polygon_1;

    // Make sure that the polygon and its copy are equal
    ASSERT_EQ(polygon_1, polygon_2);

    // Change the first node of the copy
    polygon_2[0] = {1, 1, 1};

    // Make sure that the polygon and its copy are no longer equal
    ASSERT_NE(polygon_1, polygon_2);

    // Make a polygon from 3 points
    auto polygon_3 = make_polygon(p1, p2, p3);

    // Since the 3 points are the same that forms the 3 segments when we close them, make sure that this new polygon
    // is equal to the first one
    ASSERT_EQ(polygon_1, polygon_3);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
