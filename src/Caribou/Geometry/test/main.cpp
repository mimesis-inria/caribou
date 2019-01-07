#include <gtest/gtest.h>

#include <Caribou/Traits.h>
#include <Caribou/Geometry/Node.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Interpolation/Triangle.h>
#include <Caribou/Geometry/Interpolation/Quad.h>

TEST(Geometry, Node) {
    using namespace caribou::geometry;

    auto p1 = Node(1, 2, 3);

    auto p2 = Node<3> {1, 2, 3};

    ASSERT_EQ(p1, p2);

    constexpr caribou::algebra::Vector<3> v (1, 2, 3);
    Node p3 (v);

    Node p4 = p3;

    p4[0] = p1[0];
    p4[1]= p2[1];
    p4[2] = p3[2];

    ASSERT_EQ(p3, p4);

    ASSERT_EQ(p1.x(), 1);

    p2.y() = 4;
    ASSERT_EQ(p2[1], 4);

    ASSERT_EQ(p4[2], 3);
    ASSERT_NE(p1, p2);

    Node<3> p8;

    // Make sure a point have minimal space taken in the memory
    ASSERT_EQ(
            sizeof(p8),
            (3 * sizeof(Node<3>::ValueType))
    );

    // Transformations test
    Node p9 = p1.translated({5,5,5}); // Does not change the value of p1
    ASSERT_EQ(Node(1, 2, 3), p1);
    ASSERT_EQ(Node(1+5, 2+5, 3+5), p9);

    p1.translate({5,5,5});
    ASSERT_EQ(Node(1+5, 2+5, 3+5), p1);
}

TEST(Geometry, Triangle) {
    using namespace caribou::geometry;
    using namespace caribou::geometry::interpolation;
    using namespace caribou::algebra;

    Triangle<2> triangle (Node<2> {0,0}, Node<2> {1, 0}, Node<2> {0, 1});
    std::array<float, 2> node {{1, 0}};
    ASSERT_EQ( triangle.node(1), node);

    // Shape functions
    ASSERT_EQ(Triangle3::L<0>(0, 0), 1);
    ASSERT_EQ(Triangle3::L<0>(1, 0), 0);
    ASSERT_EQ(Triangle3::L<0>(0, 1), 0);
    ASSERT_EQ(Triangle3::L<1>(1, 0), 1);
    ASSERT_EQ(Triangle3::L<2>(0, 1), 1);

    Vector <3, float> shapes {Triangle3::L<0>(0, 0), Triangle3::L<1>(0, 0), Triangle3::L<2>(0, 0)};
    ASSERT_EQ(Triangle3::N(0,0), shapes);

    // Shape function derivatives
    caribou::algebra::Vector<2, float> derivatives[3] {{-1, -1}, {1, 0}, {0, 1}};
    ASSERT_EQ(Triangle3::dL<0>(0, 0), derivatives[0]);
    ASSERT_EQ(Triangle3::dL<1>(1, 0), derivatives[1]);
    ASSERT_EQ(Triangle3::dL<2>(0, 1), derivatives[2]);

    Matrix sderivatives  {
            Triangle3::dL<0>(0, 0),
            Triangle3::dL<1>(0, 0),
            Triangle3::dL<2>(0, 0),
    };
    ASSERT_EQ(triangle.dN(0,0), sderivatives);

    // Jacobian
    Triangle<2> triangle2D (Node<2>{50, 50}, Node<2>{60, 50}, Node<2>{55, 55});
    Matrix J = triangle2D.jacobian(1/3, 1/3);
    ASSERT_EQ(J.determinant()*0.5, 25);

    // Some properties
    caribou::algebra::Vector<2, float> center {1/3., 1/3.};
    ASSERT_EQ(triangle.center(), center);

    ASSERT_EQ(triangle.area(), 0.5);

    Triangle<3> triangle1 (Node<3> {10,10,10}, Node<3> {20, 10, 10}, Node<3> {15, 20, 10});
    ASSERT_EQ(triangle1.area(), 50);

}

TEST(Geometry, Quad) {
    using namespace caribou::geometry;
    using namespace caribou::geometry::interpolation;
    using namespace caribou::algebra;

    Quad<2, interpolation::Quad4> quad (Node<2> {-1, -1}, Node<2> {+1, -1}, Node<2>{+1, +1}, Node<2>{-1, +1});

    std::array<float, 2> node {{+1, +1}};
    ASSERT_EQ( quad.node(2), node);

    // Shape functions
    ASSERT_EQ(Quad4::L<0>(-1, -1), 1);
    ASSERT_EQ(Quad4::L<0>(-1, +1), 0);
    ASSERT_EQ(Quad4::L<0>(+1, -1), 0);
    ASSERT_EQ(Quad4::L<0>(+1, +1), 0);
    ASSERT_EQ(Quad4::L<1>(+1, -1), 1);
    ASSERT_EQ(Quad4::L<2>(+1, +1), 1);
    ASSERT_EQ(Quad4::L<3>(-1, +1), 1);

    // Shape function derivatives
    caribou::algebra::Vector<2, float> derivatives[4] {{-0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {-0.5, 0.5}};
    ASSERT_EQ(Quad4::dL<0>(-1, -1), derivatives[0]);
    ASSERT_EQ(Quad4::dL<1>(+1, -1), derivatives[1]);
    ASSERT_EQ(Quad4::dL<2>(+1, +1), derivatives[2]);
    ASSERT_EQ(Quad4::dL<3>(-1, +1), derivatives[3]);

    // Jacobian
    Quad<2> quad2D (Node<2>{50, 50}, Node<2>{55, 50}, Node<2>{55, 55}, Node<2>{50, 55});
    Matrix J = quad2D.jacobian(-sqrt(3)/3, sqrt(3)/3);
    ASSERT_EQ(J.determinant()*4, 25);


    // Interpolation

    auto f = [] (const auto & u, const auto & v) -> float {
        return 5 + 2*u + 3*v;
    };

    ASSERT_EQ(f(0,0), quad.interpolate_at_local_position({0, 0}, f(-1, -1), f(1, -1), f(1, 1), f(-1, 1)));

}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
