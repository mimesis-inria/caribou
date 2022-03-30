#pragma once

#include <Eigen/Dense>
#include <Caribou/macros.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/RectangularHexahedron20.h>

TEST(Hexahedron, RectangularQuadratic) {
    using namespace caribou;
    using Hexahedron = caribou::geometry::Hexahedron20;
    using RectangularHexahedron = caribou::geometry::RectangularHexahedron20;
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
        Eigen::Matrix<FLOATING_POINT_TYPE, 8, 3> transformed_nodes;
        for (unsigned int node_id = 0; node_id < 8; ++node_id) {
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

}

