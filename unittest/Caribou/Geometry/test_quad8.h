#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <array>
#include <Caribou/macros.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Quad8.h>

TEST(Quad, Quadratic) {
    using namespace caribou;
    using std::cos;
    using std::sin;
    const double pi = std::acos(-1);

    // 2D
    {
        using Rotation = Eigen::Matrix<FLOATING_POINT_TYPE, 2, 2>;
        using Triangle = caribou::geometry::Triangle6<_2D>;
        using Quad = caribou::geometry::Quad8<_2D>;

        using LocalCoordinates = Quad::LocalCoordinates;
        using WorldCoordinates = Quad::WorldCoordinates;

        Quad q (
            WorldCoordinates(-5, -5), WorldCoordinates(+5, -5),
            WorldCoordinates(+10, +5), WorldCoordinates(-10, +5)
        );

        EXPECT_EQ(q.number_of_boundary_elements(), static_cast<UNSIGNED_INTEGER_TYPE>(4));

        // Center
        EXPECT_DOUBLE_EQ(q.center()[0], 0);
        EXPECT_DOUBLE_EQ(q.center()[1], 0);

        // Inverse transformation
        for (const auto & gauss_node : q.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, q.local_coordinates(q.world_coordinates(gauss_node.position)), 1e-5);
        }
        for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < q.number_of_nodes();++node_id) {
            EXPECT_MATRIX_NEAR(q.node(node_id), q.world_coordinates(q.local_coordinates(q.node(node_id))), 1e-5);
        }

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                q.local_coordinates(q.center()),
                q.local_coordinates(q.world_coordinates({-1., 0.5})),
                q.local_coordinates(q.world_coordinates({+1, -0.5})),
                q.local_coordinates(q.world_coordinates({0.5, +1})),
                q.local_coordinates(q.world_coordinates({-0.5, -1}))
            };
            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < q.number_of_nodes();++node_id) {
                inside_points.push_back(q.local_coordinates(q.node(node_id)));
            }
            for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < q.number_of_gauss_nodes();++gauss_node_id) {
                inside_points.push_back(q.gauss_node(gauss_node_id).position);
            }
            for (const auto & p : inside_points) {
                ASSERT_TRUE(q.contains_local(p, epsilon)) <<
                                                          "Local point [" << p[0] << ", " << p[1] << "] is found outside the element, but it should be inside.";
            }
            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {
                q.local_coordinates(q.node(0)) - LocalCoordinates(epsilon*1.1, 0),
                q.local_coordinates(q.node(1)) + LocalCoordinates(epsilon*1.1, 0),
                q.local_coordinates(q.node(2)) + LocalCoordinates(epsilon*1.1, epsilon*1.1),
                q.local_coordinates(q.node(3)) + LocalCoordinates(0, epsilon*1.1),
            };
            for (const auto & p : outside_points) {
                ASSERT_FALSE(q.contains_local(p, epsilon)) <<
                                                           "Local point [" << p[0] << ", " << p[1] << "] is found inside the element, but it should be outside.";
            }
        }

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 8, 1> values;
        values << p2(q.node(0)), p2(q.node(1)), p2(q.node(2)), p2(q.node(3)),
                  p2(q.node(4)), p2(q.node(5)), p2(q.node(6)), p2(q.node(7));
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_NEAR(q.interpolate(x, values), p2(q.world_coordinates(x)), 1e-10);
        }

        // Integration
        // Quad
        FLOATING_POINT_TYPE numerical_solution_quad = 0;
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = q.jacobian(x).determinant();
            numerical_solution_quad += p2(q.world_coordinates(x)) * w * detJ;
        }

        // Two triangles
        FLOATING_POINT_TYPE numerical_solution_triangles = 0;
        Triangle t1 (q.node(0), q.node(1), q.node(3));
        Triangle t2 (q.node(1), q.node(2), q.node(3));
        for (const auto & t : {t1, t2}) {
            for (const auto & gauss_node : t.gauss_nodes()) {
                const auto x = gauss_node.position;
                const auto w = gauss_node.weight;
                const auto detJ = t.jacobian(x).determinant();
                numerical_solution_triangles += p2(t.world_coordinates(x)) * w * detJ;
            }
        }

        EXPECT_DOUBLE_EQ(numerical_solution_quad, numerical_solution_triangles);

        // Frame extraction
        Quad initial_quad (WorldCoordinates(-25, -25), WorldCoordinates( 25, -25), WorldCoordinates(25, 25), WorldCoordinates(-25, 25));

        bool test_failed = false;
        std::array<FLOATING_POINT_TYPE, 8> angles {pi/6, pi/4, pi/3, pi/2, 3*pi/4, pi, 5*pi/4, 7*pi/4};
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 2*angles.size()> results; // Should be 8 2x2 Identity matrices
        unsigned int angle_id = 0;
        for (const auto & a : angles) {
            Rotation R;
            R << cos(a), -sin(a),
                sin(a),  cos(a);

            Eigen::Matrix<FLOATING_POINT_TYPE, 8, 2> rotated_nodes;
            for (unsigned int i = 0; i < 8; ++i)
                rotated_nodes.row(i) = R*initial_quad.node(i);

            Quad rotated_quad(rotated_nodes);
            Rotation F = rotated_quad.frame({0, 0});
            results.block<2,2>(0, 2*angle_id) = F.transpose()*R;
            Rotation diff = R - F;
            if (not (IN_CLOSED_INTERVAL(-1e-15, diff.minCoeff(), 1e-15) and
                     IN_CLOSED_INTERVAL(-1e-15, diff.maxCoeff(), 1e-15))) {
                test_failed = true;
            }
            ++angle_id;
        }
        Eigen::IOFormat clean(2, 0, ", ", "\n", "[", "]");
        EXPECT_FALSE(test_failed) << results.format(clean);
    }

    // 3D
    {
        using Triangle = caribou::geometry::Triangle6<_3D>;
        using Quad = caribou::geometry::Quad8<_3D>;

        using LocalCoordinates = Quad::LocalCoordinates;
        using WorldCoordinates = Quad::WorldCoordinates;

        Quad q (
            WorldCoordinates(-5, -53./15, -53./15), WorldCoordinates(+5, -53./15, -53./15),
            WorldCoordinates(+10,+53./15, +53./15), WorldCoordinates(-10,+53./15, +53./15)
        );

        EXPECT_EQ(q.number_of_boundary_elements(), static_cast<UNSIGNED_INTEGER_TYPE>(4));

        // Center
        EXPECT_DOUBLE_EQ(q.center()[0], 0);
        EXPECT_DOUBLE_EQ(q.center()[1], 0);
        EXPECT_DOUBLE_EQ(q.center()[2], 0);

        // Inverse transformation
        for (const auto & gauss_node : q.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, q.local_coordinates(q.world_coordinates(gauss_node.position)), 1e-5);
        }
        for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < q.number_of_nodes();++node_id) {
            EXPECT_MATRIX_NEAR(q.node(node_id), q.world_coordinates(q.local_coordinates(q.node(node_id))), 1e-5);
        }

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                q.local_coordinates(q.center()),
                q.local_coordinates(q.world_coordinates({-1., 0.5})),
                q.local_coordinates(q.world_coordinates({+1, -0.5})),
                q.local_coordinates(q.world_coordinates({0.5, +1})),
                q.local_coordinates(q.world_coordinates({-0.5, -1}))
            };
            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < q.number_of_nodes();++node_id) {
                inside_points.push_back(q.local_coordinates(q.node(node_id)));
            }
            for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < q.number_of_gauss_nodes();++gauss_node_id) {
                inside_points.push_back(q.gauss_node(gauss_node_id).position);
            }
            for (const auto & p : inside_points) {
                ASSERT_TRUE(q.contains_local(p, epsilon)) <<
                                                          "Local point [" << p[0] << ", " << p[1] << "] is found outside the element, but it should be inside.";
            }
            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {
                q.local_coordinates(q.node(0)) - LocalCoordinates(epsilon*1.1, 0),
                q.local_coordinates(q.node(1)) + LocalCoordinates(epsilon*1.1, 0),
                q.local_coordinates(q.node(2)) + LocalCoordinates(epsilon*1.1, epsilon*1.1),
                q.local_coordinates(q.node(3)) + LocalCoordinates(0, epsilon*1.1),
            };
            for (const auto & p : outside_points) {
                ASSERT_FALSE(q.contains_local(p, epsilon)) <<
                                                           "Local point [" << p[0] << ", " << p[1] << "] is found inside the element, but it should be outside.";
            }
        }

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 8, 1> values;
        values << p2(q.node(0)), p2(q.node(1)), p2(q.node(2)), p2(q.node(3)),
                  p2(q.node(4)), p2(q.node(5)), p2(q.node(6)), p2(q.node(7));
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(q.interpolate(x, values), p2(q.world_coordinates(x)));
        }

        // Integration
        // Quad
        FLOATING_POINT_TYPE numerical_solution_quad = 0;
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J = q.jacobian(x);
            const auto detJ = J.col(0).cross(J.col(1)).norm();
            numerical_solution_quad += p2(q.world_coordinates(x)) * w * detJ;
        }

        // Two triangles
        FLOATING_POINT_TYPE numerical_solution_triangles = 0;
        Triangle t1 (q.node(0), q.node(1), q.node(3));
        Triangle t2 (q.node(1), q.node(2), q.node(3));
        for (const auto & t : {t1, t2}) {
            for (const auto & gauss_node : t.gauss_nodes()) {
                const auto x = gauss_node.position;
                const auto w = gauss_node.weight;
                const auto J = t.jacobian(x);
                const auto detJ = J.col(0).cross(J.col(1)).norm();
                numerical_solution_triangles += p2(t.world_coordinates(x)) * w * detJ;
            }
        }

        EXPECT_DOUBLE_EQ(numerical_solution_quad, numerical_solution_triangles);
    }
}
