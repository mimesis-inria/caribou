#pragma once

#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Triangle.h>

TEST(Triangle, Linear) {
    using namespace caribou;

    // 2D
    {
        using Triangle = caribou::geometry::Triangle<_2D, Linear>;

        using LocalCoordinates = Triangle::LocalCoordinates;
        using WordCoordinates = Triangle::WorldCoordinates;

        Triangle t (
            WordCoordinates({50, 50}),
            WordCoordinates({60, 50}),
            WordCoordinates({55, 55})
        );

        EXPECT_EQ(t.number_of_boundary_elements(), static_cast<UNSIGNED_INTEGER_TYPE>(3));

        auto J = t.jacobian(LocalCoordinates {1/3., 1/3.} );
        EXPECT_DOUBLE_EQ(J.determinant()*0.5, 25);

        // Center
        WordCoordinates center = (t.node(0) + t.node(1) + t.node(2)) / 3.;

        EXPECT_DOUBLE_EQ(t.center()[0], center[0]);
        EXPECT_DOUBLE_EQ(t.center()[1], center[1]);

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                t.local_coordinates(t.node(0)),
                t.local_coordinates(t.node(1)),
                t.local_coordinates(t.node(2)),
                t.local_coordinates(t.center()),
                t.gauss_node(0).position,
                t.local_coordinates(t.world_coordinates({-epsilon/2., 0.5})),
                t.local_coordinates(t.world_coordinates({0.5, -epsilon/2}))
            };
            for (const auto & p : inside_points) {
                ASSERT_TRUE(t.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << "] is found outside the triangle, but it should be inside.";
            }
            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {
                t.local_coordinates(t.node(0)) - LocalCoordinates(epsilon*1.1, 0),
                t.local_coordinates(t.node(1)) + LocalCoordinates(epsilon*1.1, 0),
                t.local_coordinates(t.node(2)) + LocalCoordinates(0, epsilon*1.1),
                t.local_coordinates(t.world_coordinates({-epsilon*1.1, 0.5})),
                t.local_coordinates(t.world_coordinates({0.5, -epsilon*1.1}))
            };
            for (const auto & p : outside_points) {
                ASSERT_FALSE(t.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << "] is found inside the triangle, but it should be outside.";
            }
        }

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p1(t.node(0)), p1(t.node(1)), p1(t.node(2)));
        for (const auto & gauss_node : t.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(t.interpolate(x, values), p1(t.world_coordinates(x)));
        }

        // Inverse transformation
        for (const auto & gauss_node : t.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, t.local_coordinates(t.world_coordinates(gauss_node.position)), 1e-5);
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : t.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = t.jacobian(x).determinant();
            numerical_solution += p1(t.world_coordinates(x)) * w * detJ;
        }
        FLOATING_POINT_TYPE a = t.node(0)[0];
        FLOATING_POINT_TYPE b = t.node(1)[0];
        FLOATING_POINT_TYPE c = t.node(2)[0];
        FLOATING_POINT_TYPE d = t.node(0)[1];
        FLOATING_POINT_TYPE e = t.node(1)[1];
        FLOATING_POINT_TYPE f = t.node(2)[1];
        auto analytic_solution = static_cast<FLOATING_POINT_TYPE>(
            ((a - c)*(d - f)*(2*a + 4*c + 3*(5 + 2*d + f)))/6 -
            ((b - c)*(6*(5 + b + c)*d + 9*d*d - 3*f*f - e*(15 + 4*b + 2*c + 3*e) - f*(15 + 2*b + 4*c + 3*e)))/6);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);
    }

    // 3D
    {
        using Triangle = caribou::geometry::Triangle<_3D, Linear>;

        using LocalCoordinates = Triangle::LocalCoordinates;
        using WordCoordinates = Triangle::WorldCoordinates;

        Triangle t (
            WordCoordinates({50, 50, 0}),
            WordCoordinates({60, 50, 0}),
            WordCoordinates({55, 55, 0})
        );

        // Center
        WordCoordinates center = (t.node(0) + t.node(1) + t.node(2)) / 3.;

        EXPECT_DOUBLE_EQ(t.center()[0], center[0]);
        EXPECT_DOUBLE_EQ(t.center()[1], center[1]);
        EXPECT_DOUBLE_EQ(t.center()[2], center[2]);

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                t.local_coordinates(t.node(0)),
                t.local_coordinates(t.node(1)),
                t.local_coordinates(t.node(2)),
                t.local_coordinates(t.center()),
                t.gauss_node(0).position,
                t.local_coordinates(t.world_coordinates({-epsilon/2., 0.5})),
                t.local_coordinates(t.world_coordinates({0.5, -epsilon/2}))
            };
            for (const auto & p : inside_points) {
                ASSERT_TRUE(t.contains_local(p, epsilon)) <<
                                                          "Local point [" << p[0] << ", " << p[1] << "] is found outside the triangle, but it should be inside.";
            }
            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {
                t.local_coordinates(t.node(0)) - LocalCoordinates(epsilon*1.1, 0),
                t.local_coordinates(t.node(1)) + LocalCoordinates(epsilon*1.1, 0),
                t.local_coordinates(t.node(2)) + LocalCoordinates(0, epsilon*1.1),
                t.local_coordinates(t.world_coordinates({-epsilon*1.1, 0.5})),
                t.local_coordinates(t.world_coordinates({0.5, -epsilon*1.1}))
            };
            for (const auto & p : outside_points) {
                ASSERT_FALSE(t.contains_local(p, epsilon)) <<
                                                           "Local point [" << p[0] << ", " << p[1] << "] is found inside the triangle, but it should be outside.";
            }
        }

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p1(t.node(0)), p1(t.node(1)), p1(t.node(2)));
        for (const auto & gauss_node : t.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(t.interpolate(x, values), p1(t.world_coordinates(x)));
        }

        // Inverse transformation
        for (const auto & gauss_node : t.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, t.local_coordinates(t.world_coordinates(gauss_node.position)), 1e-5);
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : t.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J = t.jacobian(x);
            const auto detJ = J.col(0).cross(J.col(1)).norm();
            numerical_solution += p1(t.world_coordinates(x)) * w * detJ;
        }
        FLOATING_POINT_TYPE a = t.node(0)[0];
        FLOATING_POINT_TYPE b = t.node(1)[0];
        FLOATING_POINT_TYPE c = t.node(2)[0];
        FLOATING_POINT_TYPE d = t.node(0)[1];
        FLOATING_POINT_TYPE e = t.node(1)[1];
        FLOATING_POINT_TYPE f = t.node(2)[1];
        auto analytic_solution = static_cast<FLOATING_POINT_TYPE>(
            ((a - c)*(d - f)*(2*a + 4*c + 3*(5 + 2*d + f)))/6 -
            ((b - c)*(6*(5 + b + c)*d + 9*d*d - 3*f*f - e*(15 + 4*b + 2*c + 3*e) - f*(15 + 2*b + 4*c + 3*e)))/6);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);
    }
}

TEST(Triangle, Quadratic) {
    using namespace caribou;
    using Triangle = caribou::geometry::Triangle<_2D, Quadratic>;


    // 2D
    {
        using LocalCoordinates = Triangle::LocalCoordinates;
        using WordCoordinates = Triangle::WorldCoordinates;

        // Jacobian
        Triangle t (
            WordCoordinates({50, 50}),
            WordCoordinates({60, 50}),
            WordCoordinates({55, 55})
        );

        auto J = t.jacobian(LocalCoordinates {1/3., 1/3.} );
        EXPECT_NEAR(J.determinant(), 50., 1e-4);

        // Center
        WordCoordinates center = (t.node(0) + t.node(1) + t.node(2)) / 3.;

        EXPECT_DOUBLE_EQ(t.center()[0], center[0]);
        EXPECT_DOUBLE_EQ(t.center()[1], center[1]);

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                t.local_coordinates(t.node(0)),
                t.local_coordinates(t.node(1)),
                t.local_coordinates(t.node(2)),
                t.local_coordinates(t.node(3)),
                t.local_coordinates(t.node(4)),
                t.local_coordinates(t.node(5)),
                t.local_coordinates(t.center()),
                t.gauss_node(0).position,
                t.local_coordinates(t.world_coordinates({-epsilon/2., 0.5})),
                t.local_coordinates(t.world_coordinates({0.5, -epsilon/2}))
            };
            for (const auto & p : inside_points) {
                ASSERT_TRUE(t.contains_local(p, epsilon)) <<
                                                          "Local point [" << p[0] << ", " << p[1] << "] is found outside the triangle, but it should be inside.";
            }
            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {
                t.local_coordinates(t.node(0)) - LocalCoordinates(epsilon*1.1, 0),
                t.local_coordinates(t.node(1)) + LocalCoordinates(epsilon*1.1, 0),
                t.local_coordinates(t.node(2)) + LocalCoordinates(0, epsilon*1.1),
                t.local_coordinates(t.world_coordinates({-epsilon*1.1, 0.5})),
                t.local_coordinates(t.world_coordinates({0.5, -epsilon*1.1}))
            };
            for (const auto & p : outside_points) {
                ASSERT_FALSE(t.contains_local(p, epsilon)) <<
                                                           "Local point [" << p[0] << ", " << p[1] << "] is found inside the triangle, but it should be outside.";
            }
        }

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 6, 1> values;
        values << p2(t.node(0)), p2(t.node(1)), p2(t.node(2)),
                  p2(t.node(3)), p2(t.node(4)), p2(t.node(5));
        for (const auto & gauss_node : t.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(t.interpolate(x, values), p2(t.world_coordinates(x)));
        }

        // Inverse transformation
        for (const auto & gauss_node : t.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, t.local_coordinates(t.world_coordinates(gauss_node.position)), 1e-5);
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : t.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = t.jacobian(x).determinant();
            numerical_solution += p1(t.world_coordinates(x)) * w * detJ;
        }
        FLOATING_POINT_TYPE a = t.node(0)[0];
        FLOATING_POINT_TYPE b = t.node(1)[0];
        FLOATING_POINT_TYPE c = t.node(2)[0];
        FLOATING_POINT_TYPE d = t.node(0)[1];
        FLOATING_POINT_TYPE e = t.node(1)[1];
        FLOATING_POINT_TYPE f = t.node(2)[1];
        auto analytic_solution = static_cast<FLOATING_POINT_TYPE>(
            ((a - c)*(d - f)*(2*a + 4*c + 3*(5 + 2*d + f)))/6 +
           -((b - c)*(6*(5 + b + c)*d + 9*d*d - 3*f*f - e*(15 + 4*b + 2*c + 3*e) - f*(15 + 2*b + 4*c + 3*e)))/6);
        EXPECT_NEAR(numerical_solution, analytic_solution, 1e-10);
    }
}