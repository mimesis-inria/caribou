#include <gtest/gtest.h>
#include <Caribou/Geometry/Point.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Polygon.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/LinearHexahedron.h>
#include <Caribou/Geometry/LinearRegularHexahedron.h>

TEST(Geometry, Point) {
    using namespace caribou::geometry;

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

    Point3D p8;

    // Make sure a point have minimal space taken in the memory
    ASSERT_EQ(
        sizeof(p8),
        (3 * sizeof(Point3D::ValueType))
    );
}

TEST(Geometry, Segment) {
    using namespace caribou::geometry;

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
    Point3D p5;
    Point3D p6;
    auto s4 = make_segment(p5, p6);
    ASSERT_LE(
        sizeof(s4),
        (sizeof(p5) + sizeof(p6) + 1)
    );

    // Test the segment contains point
    auto s5 = make_segment({0, 0, 0}, {1, 1, 1});
    ASSERT_TRUE(s5.contains({0, 0, 0}));
    ASSERT_TRUE(s5.contains({1, 1, 1}));
    ASSERT_TRUE(s5.contains({0.5, 0.5, 0.5}));

    ASSERT_FALSE(s5.contains({1.1, 1.1, 1.1}));
    ASSERT_FALSE(s5.contains({-0.1, -0.1, -0.1}));
}

TEST(Geometry, Polygon) {
    using namespace caribou::geometry;

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

    Point3D p1 = {0, 0, 0};
    Point3D p2 = {1, 0, 0};
    Point3D p3 = {0, 1, 0};

    // Creates a triangle from 3 points
    auto t1 = make_triangle(p1, p2, p3);

    // Creates the same triangle with 3 lists initializer
    auto t2 = make_triangle({0,0,0}, {1,0,0}, {0,1,0});

    // Creates the same triangle with a list initializer of 3 lists initializer
    auto t3 = make_triangle({{0,0,0}, {1,0,0}, {0,1,0}});

    Triangle<3> t4 = t3;

    // Makes sure they are all equals
    ASSERT_EQ(t1, t2);
    ASSERT_EQ(t1, t3);
    ASSERT_EQ(t1, t4);

    // Test the center of a triangle
    ASSERT_GT(0.0001, (t4.center() - make_point(0.333333, 0.333333, 0)).length());

    // Test the normal of a triangle
    ASSERT_EQ(t3.normal(), caribou::algebra::Vector<3>(0, 0, 1));
}

TEST(Geometry, Quad) {
    using namespace caribou::geometry;

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
    using namespace caribou::algebra;

    using Vector = LinearHexahedron::VectorType;
    using Float = LinearHexahedron::Float;

    // Constructors test
    Point3D p1, p2, p3, p4, p5, p6, p7, p8;
    LinearHexahedron hexa_1 ({p1, p2, p3, p4, p5, p6, p7, p8});

    Vector v1, v2, v3, v4, v5, v6, v7, v8;
    LinearHexahedron hexa_2 ({v1, v2, v3, v4, v5, v6, v7, v8});

    LinearRegularHexahedron hexa_3({1,1,1}, {5, 5, 5});

    // Make sure a hexahedron does not take anymore space than 8 nodes
    ASSERT_EQ(8 * sizeof(p1), sizeof(hexa_1));
    ASSERT_EQ(2 * sizeof(p1), sizeof(hexa_3));

    ASSERT_EQ(Vector({1,1,1}), hexa_3.node(0));
    ASSERT_EQ(Vector({6,1,1}), hexa_3.node(1));
    ASSERT_EQ(Vector({6,6,1}), hexa_3.node(2));
    ASSERT_EQ(Vector({1,6,1}), hexa_3.node(3));
    ASSERT_EQ(Vector({1,1,6}), hexa_3.node(4));
    ASSERT_EQ(Vector({6,1,6}), hexa_3.node(5));
    ASSERT_EQ(Vector({6,6,6}), hexa_3.node(6));
    ASSERT_EQ(Vector({1,6,6}), hexa_3.node(7));

    // Create a base hexahedron by using the base constructor without any arguments
    LinearHexahedron base_hexa;
    auto edge_0 = base_hexa.edge(0);
    ASSERT_EQ(2, edge_0.length());

    auto face_0 = base_hexa.face(0);
    ASSERT_EQ(Vector(0, 0, 1), face_0.normal());

    // Shape functions
    // World coordinates should be equal to local coordinates if the hexahedron is elemental (-1 to 1 for x, y and z axis)
    ASSERT_EQ(Vector(0, 0, 0), LinearHexahedron().from_local_coordinate(Vector(0, 0, 0)));

    // Simple geometric tests
    // Makes sure the jacobian of any position inside a linear hexahedron is the same
    ASSERT_EQ(
            (int) LinearHexahedron().Jacobian(-1./sqrt(3.0), -1./sqrt(3.0), -1./sqrt(3.0)).determinant(),
            (int) LinearHexahedron().Jacobian(+1./sqrt(3.0), +1./sqrt(3.0), +1./sqrt(3.0)).determinant()
    );

    // Makes sure the 8 times the determinant of the jacobian of a linear hexahedron is equal to the volume
    ASSERT_EQ(8*8*8, LinearHexahedron().scale(4).Jacobian(0,0,0).determinant()*8);

    // Quadrature
    Float volume = LinearHexahedron().gauss_quadrature(
            (Float) 0,
            [](const LinearHexahedron & /*h*/, const Float & /*xi*/, const Float & /*eta*/, const Float & /*zeta*/) -> Float {
                return 1;
            }
    );

    // Makes sure the integral of "1" on the hexahedron equals its volume
    ASSERT_FLOAT_EQ(8, volume);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
