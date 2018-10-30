#include <gtest/gtest.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Hexahedron.h>

TEST(Geometry, Point) {
    using namespace caribou::geometry;
    using Point3D = Point3D<>;

    auto p1 = make_point(1, 2, 3);
    auto p2 = make_point(1, 2, 3);

    ASSERT_EQ(p1, p2);

    Point3D p3({1, 2, 3});
    Point3D p4;

    p4[0] = p1[0];
    p4[1]= p2[1];
    p4[2] = p3[2];

    ASSERT_EQ(p3, p4);

    auto p5 = p4;

    ASSERT_EQ(p5.x(), 1);

    p5.y() = 4;
    ASSERT_EQ(p5[1], 4);

    ASSERT_EQ(p5[2], 3);
    ASSERT_EQ(p1, p2);

    caribou::geometry::Point3D<caribou::algebra::Vector<3, unsigned char>> p8;

    // Make sure a point have minimal space taken in the memory
    ASSERT_LE(
        sizeof(p8),
        (3 * sizeof(unsigned char) + 1)
    );
}

TEST(Geometry, Segment) {
    using namespace caribou::geometry;
    using Point3D = Point3D<>;

    Point3D p1(0, 0, 0);
    Point3D p2 = {1, 1, 1};
    auto s1 =  make_segment(p1, p2);

    Point3D p3(0, 0, 0);
    Point3D p4 = {1, 1, 1};
    auto s2 = make_segment(p3, p4);

    ASSERT_EQ(s1, s2);
    ASSERT_EQ(s2, s1);

    auto s3 = make_segment({0, 0, 0}, {1, 1, 1});
    ASSERT_EQ(s1, s3);

    // Make sure a segment have minimal space taken in the memory
    caribou::geometry::Point3D<caribou::algebra::Vector<3, unsigned char>> p5;
    caribou::geometry::Point3D<caribou::algebra::Vector<3, unsigned char>> p6;
    auto s4 = make_segment(p5, p6);
    ASSERT_LE(
        sizeof(s4),
        (sizeof(p5) + sizeof(p6) + 1)
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

    // Test polygon creation with initializer lists
    auto polygon_4 = make_polygon({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}});
    ASSERT_EQ(polygon_3, polygon_4);
}

TEST(Geometry, Triangle) {
    using namespace caribou::geometry;
    using Point3D = Point3D<>;

    Point3D p1 = {0, 0, 0};
    Point3D p2 = {1, 0, 0};
    Point3D p3 = {0, 1, 0};

    // Creates a triangle from 3 points
    auto t1 = make_triangle(p1, p2, p3);

    // Creates the same triangle with 3 lists initializer
    auto t2 = make_triangle({0,0,0}, {1,0,0}, {0,1,0});

    // Creates the same triangle with a list initializer of 3 lists initializer
    auto t3 = make_triangle({{0,0,0}, {1,0,0}, {0,1,0}});

    // Makes sure they are all equals
    ASSERT_EQ(t1, t2);
    ASSERT_EQ(t1, t3);

    ASSERT_EQ(t3.normal(), caribou::algebra::Vector<3>(0, 0, 1));
}

TEST(Geometry, Quad) {
    using namespace caribou::geometry;
    using Point3D = Point3D<>;

    Point3D p1 = {1, -1, 0};
    Point3D p2 = {-1, -1, 0};
    Point3D p3 = {1, 1, 0};
    Point3D p4 = {-1, 1, 0};

    // Creates a quad from 4 points
    auto q1 = make_quad(p1, p2, p3, p4);

    // Creates the same quad with 3 lists initializer
    auto q2 = make_quad({1, -1, 0}, {-1, -1, 0}, {1, 1, 0}, {-1, 1, 0});

    // Creates the same quad with a list initializer of 4 lists initializer
    auto q3 = make_quad({ {1, -1, 0}, {-1, -1, 0}, {1, 1, 0}, {-1, 1, 0} });

    // Makes sure they are all equals
    ASSERT_EQ(q1, q2);
    ASSERT_EQ(q1, q3);

    ASSERT_EQ(q3.normal(), caribou::algebra::Vector<3>(0, 0, -1));
}

TEST(Geometry, Hexahedron) {
    using namespace caribou::geometry;
    using Point3D = Point3D<>;

    // Constructors test
    Point3D p1, p2, p3, p4, p5, p6, p7, p8;
    Hexahedron<8, Point3D::VectorType> hexa_1 ({p1, p2, p3, p4, p5, p6, p7, p8});

    Point3D::VectorType v1, v2, v3, v4, v5, v6, v7, v8;
    RegularLinearHexahedron<Point3D::VectorType> hexa_2 ({v1, v2, v3, v4, v5, v6, v7, v8});

    // Make sure a hexahedron does not take anymore space than 8 nodes
    ASSERT_EQ(8 * sizeof(p1), sizeof(hexa_1));

    // Create a base hexahedron by using the base constructor without any arguments
    RegularLinearHexahedron<Point3D::VectorType> base_hexa;
    auto edge_0 = base_hexa.edge(0);
    ASSERT_EQ(2, edge_0.length());
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
