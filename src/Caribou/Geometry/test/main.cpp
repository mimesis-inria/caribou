#include <gtest/gtest.h>

#include <Caribou/config.h>
#include <Caribou/Traits.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Quad.h>
//#include <Caribou/Geometry/Hexahedron.h>
//#include <Caribou/Geometry/RectangularHexahedron.h>

//TEST(Geometry, Node) {
//    using namespace caribou::geometry;
//
//    auto p1 = Node(1, 2, 3);
//
//    auto p2 = Node<3> {1, 2, 3};
//
//    ASSERT_EQ(p1, p2);
//
//    constexpr caribou::algebra::Vector<3> v (1, 2, 3);
//    Node p3 (v);
//
//    Node p4 = p3;
//
//    p4[0] = p1[0];
//    p4[1]= p2[1];
//    p4[2] = p3[2];
//
//    ASSERT_EQ(p3, p4);
//
//    ASSERT_EQ(p1.x(), 1);
//
//    p2.y() = 4;
//    ASSERT_EQ(p2[1], 4);
//
//    ASSERT_EQ(p4[2], 3);
//    ASSERT_NE(p1, p2);
//
//    Node<3> p8;
//
//    // Make sure a point have minimal space taken in the memory
//    ASSERT_EQ(
//            sizeof(p8),
//            (3 * sizeof(Node<3>::ValueType))
//    );
//
//    // Transformations test
//    Node p9 = p1.translated({5,5,5}); // Does not change the value of p1
//    ASSERT_EQ(Node(1, 2, 3), p1);
//    ASSERT_EQ(Node(1+5, 2+5, 3+5), p9);
//
//    p1.translate({5,5,5});
//    ASSERT_EQ(Node(1+5, 2+5, 3+5), p1);
//}

TEST(Geometry, Segment) {
    using namespace caribou::geometry;
    using namespace caribou::geometry::interpolation;

    // Shape functions
    {
        using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE,1,1>;
        ASSERT_EQ(Segment2::L(LocalCoordinates(-1))[0], 1);
        ASSERT_EQ(Segment2::L(LocalCoordinates(-1))[1], 0);
        ASSERT_EQ(Segment2::L(LocalCoordinates(+1))[0], 0);
        ASSERT_EQ(Segment2::L(LocalCoordinates(+1))[1], 1);
    }

    // 1D
    {
        using LocalCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE,1,1>;
        using WordCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 1, 1>;
        Segment<1> segment(-5.5, 1.1);
        ASSERT_EQ(segment.node(0), WordCoordinates(-5.5));
        ASSERT_EQ(segment.node(1), WordCoordinates(1.1));

        ASSERT_FLOAT_EQ(segment.center()[0], -2.2);

        Eigen::Matrix<FLOATING_POINT_TYPE, 1, 1> jacobian(-5.5*-1. + 1.1*1);
        ASSERT_EQ(segment.jacobian(LocalCoordinates(0.)), jacobian);
    }

    // 2D
    {
        using WordCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1>;
        WordCoordinates node_0 {-1.5, -1.5};
        WordCoordinates node_1 {5.5, 5.5};
        Segment<2> segment(node_0, node_1);
        auto center_node = node_0 + (node_1 - node_0).normalized()*(node_1-node_0).norm()/2.;
        ASSERT_FLOAT_EQ(segment.center()[0], center_node[0]);
        ASSERT_FLOAT_EQ(segment.center()[1], center_node[1]);
    }

    // 3D
    {
        using WordCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>;
        WordCoordinates node_0 {-1.5, -1.5, -5.2};
        WordCoordinates node_1 { 5.5,  5.5, 54.3};
        Segment<3> segment(node_0, node_1);
        auto center_node = node_0 + (node_1 - node_0).normalized()*(node_1-node_0).norm()/2.;
        ASSERT_FLOAT_EQ(segment.center()[0], center_node[0]);
        ASSERT_FLOAT_EQ(segment.center()[1], center_node[1]);
        ASSERT_FLOAT_EQ(segment.center()[2], center_node[2]);
    }
}

TEST(Geometry, Triangle) {
    using namespace caribou::geometry;
    using namespace caribou::geometry::interpolation;

    using LocalCoordinates = Triangle3::LocalCoordinates;

    // Shape functions
    ASSERT_EQ(Triangle3::L({0, 0})[0], 1);
    ASSERT_EQ(Triangle3::L({1, 0})[0], 0);
    ASSERT_EQ(Triangle3::L({0, 1})[0], 0);
    ASSERT_EQ(Triangle3::L({1, 0})[1], 1);
    ASSERT_EQ(Triangle3::L({0, 1})[2], 1);

    // 2D
    {
        using WordCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1>;

        // Jacobian
        Triangle<2> t (
            WordCoordinates({50, 50}),
            WordCoordinates({60, 50}),
            WordCoordinates({55, 55})
            );

        auto J = t.jacobian(LocalCoordinates {1/3., 1/3.} );
        ASSERT_EQ(J.determinant()*0.5, 25);

        // Center
        WordCoordinates center = (t.node(0) + t.node(1) + t.node(2)) / 3.;

        ASSERT_FLOAT_EQ(t.center()[0], center[0]);
        ASSERT_FLOAT_EQ(t.center()[1], center[1]);

        // Area
        ASSERT_FLOAT_EQ(t.area(), 25.);
    }

    // 3D
    {
        using WordCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>;

        // Jacobian
        Triangle<3> t (
            WordCoordinates({50, 50, 33}),
            WordCoordinates({60, 50, 21}),
            WordCoordinates({55, 55, -4})
        );

        // Center
        WordCoordinates center = (t.node(0) + t.node(1) + t.node(2)) / 3.;

        ASSERT_FLOAT_EQ(t.center()[0], center[0]);
        ASSERT_FLOAT_EQ(t.center()[1], center[1]);
        ASSERT_FLOAT_EQ(t.center()[2], center[2]);

        // Area
        ASSERT_FLOAT_EQ(t.area(), 159.84367);
    }

}

TEST(Geometry, Quad) {
    using namespace caribou::geometry;
    using namespace caribou::geometry::interpolation;

    using LocalCoordinates = Triangle3::LocalCoordinates;

    // Shape functions
    ASSERT_EQ(Quad4::L({-1, -1})[0], 1);
    ASSERT_EQ(Quad4::L({-1, +1})[0], 0);
    ASSERT_EQ(Quad4::L({+1, -1})[0], 0);
    ASSERT_EQ(Quad4::L({+1, +1})[0], 0);
    ASSERT_EQ(Quad4::L({+1, -1})[1], 1);
    ASSERT_EQ(Quad4::L({+1, +1})[2], 1);
    ASSERT_EQ(Quad4::L({-1, +1})[3], 1);

    // 2D
    {
        using WordCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1>;

        // Jacobian
        Quad<2> quad2D (
            WordCoordinates({50, 50}),
            WordCoordinates({55, 50}),
            WordCoordinates({55, 55}),
            WordCoordinates({50, 55})
            );
        auto J = quad2D.jacobian(LocalCoordinates {-sqrt(3)/3, sqrt(3)/3 });
        ASSERT_EQ(J.determinant()*4, 25);
    }

//
//
//    // Interpolation
//
//    auto f = [] (const auto & u, const auto & v) -> FLOATING_POINT_TYPE {
//        return 5 + 2*u + 3*v;
//    };
//
//    ASSERT_EQ(f(0,0), quad.interpolate_at_local_position(LocalCoordinates {0, 0}, f(-1, -1), f(1, -1), f(1, 1), f(-1, 1)));

}

//TEST(Geometry, Hexahedron) {
//    using namespace caribou::geometry;
//    using namespace caribou::geometry::interpolation;
//    using namespace caribou::algebra;
//
//    Hexahedron<interpolation::Hexahedron8> hexa;
//    const auto I = Matrix<3,3>::Identity();
//
//    Node<3> node {{1, 1, -1}};
//    ASSERT_EQ(hexa.node(2), node);
//
//    double r = hexa.gauss_quadrature([](const Hexahedron<Hexahedron8> & /*h*/, const auto & /*local_coordinates*/) {
//        return 1.;
//    });
//
//    ASSERT_FLOAT_EQ(r, 8);
//
//    ASSERT_EQ(hexa.frame(), I);
//
//    // Rotate the hexa by R and extract the resulting frame
//    const Matrix R ({
//            {0.44480652434057482703, 0.49694411802952204171, 0.74511321251197670801},
//            {0.60593116938112601133, 0.44567310019856992698, -0.6589559209323614386},
//            {-0.65954118436719233465, 0.74459521306227582915, -0.10287562786328596776}
//    });
//
//    for (std::size_t i = 0; i<8; i++)
//        hexa.node(i) = R*hexa.node(i);
//
//    const auto frame = hexa.frame();
//    for (std::size_t i = 0; i < 3; ++i)
//        for (std::size_t j = 0; j < 3; ++j)
//            ASSERT_FLOAT_EQ(R(i,j), frame(i,j));
//
//}
//
//TEST(Geometry, RectangularHexahedron) {
//    using namespace caribou::geometry;
//    using namespace caribou::geometry::interpolation;
//    using namespace caribou::algebra;
//
//    {
//        RectangularHexahedron<interpolation::Hexahedron8> hexa;
//        const auto I = Matrix<3, 3>::Identity();
//
//        Node<3> node{{1, 1, -1}};
//        ASSERT_EQ(hexa.node(2), node);
//
//        double r = hexa.gauss_quadrature(
//            [](const RectangularHexahedron<Hexahedron8> & /*h*/, const auto & /*local_coordinates*/) {
//                return 1.;
//            });
//
//        ASSERT_FLOAT_EQ(r, 8);
//
//        ASSERT_EQ(hexa.frame(), I);
//    }
//
//    {
//        using Hexa = RectangularHexahedron<interpolation::Hexahedron8>;
//        using Segment = Segment<3>;
//        Hexa::NodeType center{25, 40, 2};
//        Hexa::Size dimensions{5,20, 10};
//
//        // Rotation from center point around x, y and then z axis => x: -25 deg, y: 30 deg, z: 90 deg
//        FLOATING_POINT_TYPE a = -25 * 180/M_PI;
//        Hexa::Mat33 Rx = {{
//            {1., 0., 0.},
//            {0., cos(a), -sin(a)},
//            {0, sin(a), cos(a)}
//        }};
//
//        a = 30 * 180/M_PI;
//        Hexa::Mat33 Ry = {{
//            {cos(a), 0., sin(a)},
//            {0., 1, 0.},
//            {-sin(a), 0, cos(a)}
//        }};
//
//
//        a = 90 * 180/M_PI;
//        Hexa::Mat33 Rz = {{
//            {cos(a), -sin(a), 0.},
//            {sin(a), cos(a), 0.},
//            {0, 0, 1}
//        }};
//
//        Hexa::Mat33 R = Rz*Ry*Rz;
//
//        Hexa h(center, dimensions, R);
//        Hexa base_hexa(Hexa::NodeType {0.,0.,0.}, Hexa::Size {2.,2.,2.});
//
//        for (unsigned int i = 0; i < 8; ++i) {
//            for (unsigned int j = 0; j < 3; ++j) {
//                ASSERT_FLOAT_EQ(h.T(base_hexa.node(i))[j], h.node(i)[j]);
//                ASSERT_FLOAT_EQ(h.Tinv(h.node(i))[j], base_hexa.node(i)[j]);
//            }
//        }
//
//        ASSERT_TRUE(h.intersects(Segment (Hexa::NodeType {0,0,0}, h.center())));
//        ASSERT_FALSE(h.intersects(Segment (Hexa::NodeType {0,0,0}, Hexa::NodeType {1,1,1})));
//
//        for (unsigned int i = 0; i < 8; ++i) {
//            ASSERT_TRUE(h.intersects(Segment(h.node(i), h.center()))) << "Segment from node " << i << " and the center of the hexa should intersect the cube.";
//            ASSERT_TRUE(h.intersects(Segment(h.center(), h.node(i)))) << "Segment from  the center of the hexa and node " << i << " should intersect the cube.";
//        }
//
//        for (unsigned int i = 0; i < 8; ++i) {
//            for (unsigned int j = 0; j < 8; ++j) {
//                if (i == j) continue;
//                ASSERT_TRUE(h.intersects(Segment(h.node(i), h.node(j)))) << "Segment from nodes " << i << " and " << j << " should intersect the cube.";
//            }
//        }
//
//    }
//
//}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
