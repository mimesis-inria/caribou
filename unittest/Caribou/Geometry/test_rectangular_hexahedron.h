#pragma once

#include <Eigen/Dense>
#include <Caribou/macros.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/RectangularHexahedron.h>

TEST(Hexahedron, RectangularLinear) {
    using namespace caribou;
    using Hexahedron = caribou::geometry::Hexahedron;
    using RectangularHexahedron = caribou::geometry::RectangularHexahedron;
    using WorldCoordinates = Hexahedron::WorldCoordinates;
    using Rotation = RectangularHexahedron::Rotation;
    using Size = RectangularHexahedron::Size;

    using std::cos;
    using std::sin;
    const double pi = std::acos(-1);

    Size H(25, 30, 75);
    WorldCoordinates center(75.4, 12.2, -45);
    FLOATING_POINT_TYPE a = 3*pi/4;
    Rotation R, Rx, Ry, Rz;
    Rx << 1,       0,       0,
          0,       cos(a), -sin(a),
          0,       sin(a),  cos(a);
    Ry << cos(a),  0,       sin(a),
          0,       1,       0,
          -sin(a), 0,       cos(a);
    Rz << cos(a),  -sin(a), 0,
          sin(a),  cos(a),  0,
          0      , 0,       1;
    R = Rz*Ry*Rx;

    // Rectangular hexahedron from a generic hexahedron
    {
        Hexahedron generic_hexa;
        Eigen::Matrix<FLOATING_POINT_TYPE, geometry::traits<Hexahedron>::NumberOfNodesAtCompileTime, 3> transformed_nodes;
        for (unsigned int node_id = 0; node_id < generic_hexa.number_of_nodes(); ++node_id) {
            const auto x = generic_hexa.node(node_id);
            transformed_nodes.row(node_id) = center + R * (x.cwiseProduct(H / 2.));
        }

        RectangularHexahedron q(transformed_nodes);
        EXPECT_MATRIX_NEAR(q.center(), center, 1e-10);
        EXPECT_MATRIX_NEAR(q.rotation(), R, 1e-10);
        EXPECT_MATRIX_NEAR(q.size(), H, 1e-10);
    }

    // Generic hexahedron from rectangular hexahedron
    {
        RectangularHexahedron rectangular_hexa(center, H, R);
        const auto transformed_nodes = rectangular_hexa.nodes();
        Hexahedron q(transformed_nodes);

        EXPECT_MATRIX_NEAR(q.center(), center, 1e-10);
        EXPECT_MATRIX_NEAR(q.frame({0, 0, 0}), R, 1e-10);
    }

    // World coordinates to local coordinates
    {
        RectangularHexahedron h(center, H, R);
        for (const auto & gauss_node : h.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_MATRIX_NEAR(x, h.local_coordinates(h.world_coordinates(x)), 1e-10);
        }
    }

    // Segment intersection
    {
        using Segment = caribou::geometry::Segment<_3D>;
        RectangularHexahedron h(center, H, R);
        EXPECT_TRUE(h.intersects(Segment (WorldCoordinates {0,0,0}, h.center())));
        EXPECT_FALSE(h.intersects(Segment (WorldCoordinates {0,0,0}, WorldCoordinates {1,1,1})));

        for (unsigned int i = 0; i < h.number_of_nodes(); ++i) {
            EXPECT_TRUE(h.intersects(Segment(h.node(i), h.center()))) << "Segment from node " << i << " and the center of the hexa should intersect the cube.";
            EXPECT_TRUE(h.intersects(Segment(h.center(), h.node(i)))) << "Segment from  the center of the hexa and node " << i << " should intersect the cube.";
        }

        for (unsigned int i = 0; i < h.number_of_nodes(); ++i) {
            for (unsigned int j = 0; j < h.number_of_nodes(); ++j) {
                EXPECT_TRUE(h.intersects(Segment(h.node(i), h.node(j)), 1e-10)) << "Segment from node " << i << " to node " << j << " should intersect the cube.";
            }
        }
    }

    // Triangle intersection
    {
        using Triangle = caribou::geometry::Triangle<_3D>;
        RectangularHexahedron h(WorldCoordinates {5.,5.,5.});
        // Inside triangle
        Triangle t1(WorldCoordinates{4,4,4}, WorldCoordinates {5,4,4}, WorldCoordinates {4.5,5,5});
        ASSERT_TRUE(h.intersects(t1));

        // Outside triangle
        Triangle t2(WorldCoordinates{4,4,8}, WorldCoordinates {5,4,8}, WorldCoordinates {4.5,5,9});
        ASSERT_FALSE(h.intersects(t2));

        // Lying on face triangle
        Triangle t3(WorldCoordinates{3,3,3}, WorldCoordinates {3,3,7}, WorldCoordinates {3,5,7});
        ASSERT_FALSE(h.intersects(t3));

        // Cut triangle
        Triangle t4(WorldCoordinates{2,2,2}, WorldCoordinates {8,2,2}, WorldCoordinates {5,8,8});
        ASSERT_TRUE(h.intersects(t4));
    }

}
