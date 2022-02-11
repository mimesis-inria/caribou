#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <array>
#include <Caribou/macros.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/RectangularQuad8.h>

TEST(Quad, RectangularQuadratic) {
    using namespace caribou;
    using std::cos;
    using std::sin;
    const double pi = std::acos(-1);

    // 2D
    {
        using Quad = caribou::geometry::Quad8<_2D>;
        using RectangularQuad = caribou::geometry::RectangularQuad8<_2D>;
        using Rotation = RectangularQuad::Rotation;
        using Size = RectangularQuad::Size;
        using WorldCoordinates = RectangularQuad::WorldCoordinates;

        Size H(25, 30);
        WorldCoordinates center(75.4, 12.2);
        Rotation R;
        R << cos(3*pi/4), -sin(3*pi/4),
        sin(3*pi/4),  cos(3*pi/4);

        // Rectangular quad from a generic quad
        {
            Quad generic_quad;
            Eigen::Matrix<FLOATING_POINT_TYPE, geometry::traits<Quad>::NumberOfNodesAtCompileTime, 2> transformed_nodes;
            for (unsigned int node_id = 0; node_id < generic_quad.number_of_nodes(); ++node_id) {
                const auto x = generic_quad.node(node_id);
                transformed_nodes.row(node_id) = center + R * (x.cwiseProduct(H / 2.));
            }

            RectangularQuad q(transformed_nodes);
            EXPECT_MATRIX_NEAR(q.center(), center, 1e-10);
            EXPECT_MATRIX_NEAR(q.rotation(), R, 1e-10);
            EXPECT_MATRIX_NEAR(q.size(), H, 1e-10);
        }

        // Generic quad from rectangular quad
        {
            RectangularQuad rectangular_quad(center, H, R);
            const auto transformed_nodes = rectangular_quad.nodes();
            Quad q(transformed_nodes);

            EXPECT_MATRIX_NEAR(q.center(), center, 1e-10);
            EXPECT_MATRIX_NEAR(q.frame({0, 0}), R, 1e-10);
        }
    }

    // 3D
    {
        using Quad = caribou::geometry::Quad8<_3D>;
        using RectangularQuad = caribou::geometry::RectangularQuad8<_3D>;
        using Rotation = RectangularQuad::Rotation;
        using WorldCoordinates = RectangularQuad::WorldCoordinates;

        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> H(25, 30, 0);
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

        // Rectangular quad from a generic quad
        {
            Quad generic_quad;
            Eigen::Matrix<FLOATING_POINT_TYPE, geometry::traits<Quad>::NumberOfNodesAtCompileTime, 3> transformed_nodes;
            for (unsigned int node_id = 0; node_id < generic_quad.number_of_nodes(); ++node_id) {
                const auto x = generic_quad.node(node_id);
                transformed_nodes.row(node_id) = center + R * (x.cwiseProduct(H / 2.));
            }

            RectangularQuad q(transformed_nodes);
            EXPECT_MATRIX_NEAR(q.center(), center, 1e-10);
            EXPECT_MATRIX_NEAR(q.rotation(), R, 1e-10);
            EXPECT_MATRIX_NEAR(q.size(), RectangularQuad::Size(H.block<2,1>(0,0)), 1e-10);
        }

        // Generic quad from rectangular quad
        {
            RectangularQuad::Size S = H.block<2,1>(0,0);
            RectangularQuad rectangular_quad(center, S, R);
            const auto transformed_nodes = rectangular_quad.nodes();
            Quad q(transformed_nodes);

            EXPECT_MATRIX_NEAR(q.center(), center, 1e-10);
            EXPECT_MATRIX_NEAR(q.frame({0, 0}), R, 1e-10);
        }
    }
}