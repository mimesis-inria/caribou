#pragma once

#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Segment3.h>

TEST(Segment, Linear) {
    using namespace caribou;

    // Shape functions
    {
        using Segment = caribou::geometry::Segment<_3D>;
        using LocalCoordinates = Segment::LocalCoordinates ;
        Segment s;
        EXPECT_EQ(s.L(LocalCoordinates(-1))[0], 1);
        EXPECT_EQ(s.L(LocalCoordinates(-1))[1], 0);
        EXPECT_EQ(s.L(LocalCoordinates(+1))[0], 0);
        EXPECT_EQ(s.L(LocalCoordinates(+1))[1], 1);
    }

// 1D
    {
        using Segment = caribou::geometry::Segment<_1D>;
        using LocalCoordinates  = Segment::LocalCoordinates;
        using WorldCoordinates  = Segment::WorldCoordinates;

        Segment segment(WorldCoordinates(-5.5), WorldCoordinates(1.1));
        EXPECT_MATRIX_NEAR(segment.node(0), WorldCoordinates(-5.5), 1e-5);
        EXPECT_MATRIX_NEAR(segment.node(1), WorldCoordinates(1.1), 1e-5);

        // Center
        EXPECT_DOUBLE_EQ(segment.center()[0], -2.2);

        // Inverse transformation
        for (const auto & gauss_node : segment.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, segment.local_coordinates(segment.world_coordinates(gauss_node.position)), 1e-5);
        }
        for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < segment.number_of_nodes();++node_id) {
            EXPECT_MATRIX_NEAR(segment.node(node_id), segment.world_coordinates(segment.local_coordinates(segment.node(node_id))), 1e-5);
        }

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                segment.local_coordinates(segment.center()),
                segment.local_coordinates(segment.world_coordinates(LocalCoordinates(-1))),
                segment.local_coordinates(segment.world_coordinates(LocalCoordinates(+1))),
            };
            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < segment.number_of_nodes();++node_id) {
                inside_points.push_back(segment.local_coordinates(segment.node(node_id)));
            }
            for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < segment.number_of_gauss_nodes();++gauss_node_id) {
                inside_points.push_back(segment.gauss_node(gauss_node_id).position);
            }
            for (const auto & p : inside_points) {
                ASSERT_TRUE(segment.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << "] is found outside the element, but it should be inside.";
            }
            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {
                segment.local_coordinates(segment.node(0)) - LocalCoordinates(epsilon*1.1),
                segment.local_coordinates(segment.node(1)) + LocalCoordinates(epsilon*1.1),
            };
            for (const auto & p : outside_points) {
                ASSERT_FALSE(segment.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << "] is found inside the element, but it should be outside.";
            }
        }

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> values (p1(segment.node(0)), p1(segment.node(1)));
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(segment.interpolate(x, values), p1(segment.world_coordinates(x)));
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = segment.jacobian(x)[0];
            numerical_solution += p1(segment.world_coordinates(x)) * w * detJ;
        }
        FLOATING_POINT_TYPE analytic_solution = (5*1.1 + (1.1*1.1)) - (-5.5*5 + (-5.5)*(-5.5));
        EXPECT_NEAR(numerical_solution, analytic_solution, 1e-10);
    }

// 2D
    {
        using Segment = caribou::geometry::Segment<_2D>;
        using WorldCoordinates  = Segment::WorldCoordinates ;

        WorldCoordinates node_0 {-1.5, -1.5};
        WorldCoordinates node_1 {5.5, 5.5};

        Segment segment(node_0, node_1);

        // Center
        WorldCoordinates center_node = node_0 + (node_1 - node_0).normalized()*(node_1-node_0).norm()/2.;
        EXPECT_DOUBLE_EQ(segment.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(segment.center()[1], center_node[1]);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> values (p1(segment.node(0)), p1(segment.node(1)));
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(segment.interpolate(x, values), p1(segment.world_coordinates(x)));
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  segment.jacobian(x);
            const auto detJ = J.norm();
            numerical_solution += p1(segment.world_coordinates(x)) * w * detJ;
        }

        const auto x0 = node_0[0];
        const auto y0 = node_0[1];
        const auto x1 = node_1[0];
        const auto y1 = node_1[1];
        const auto d = (node_1 - node_0).norm();
        FLOATING_POINT_TYPE analytic_solution = d * (x0 + x1 + 3*y0/2. + 3*y1/2. + 5);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);
    }

// 3D
    {
        using Segment = caribou::geometry::Segment<_3D>;
        using WorldCoordinates  = Segment::WorldCoordinates ;

        WorldCoordinates node_0 {-1.5, -1.5, -5.2};
        WorldCoordinates node_1 { 5.5,  5.5, 54.3};

        Segment segment(node_0, node_1);

        // Center
        WorldCoordinates center_node = node_0 + (node_1 - node_0).normalized()*(node_1-node_0).norm()/2.;
        EXPECT_DOUBLE_EQ(segment.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(segment.center()[1], center_node[1]);
        EXPECT_DOUBLE_EQ(segment.center()[2], center_node[2]);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 1> values (p1(segment.node(0)), p1(segment.node(1)));
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(segment.interpolate(x, values), p1(segment.world_coordinates(x)));
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  segment.jacobian(x);
            const auto detJ = J.norm();
            numerical_solution += p1(segment.world_coordinates(x)) * w * detJ;
        }

        const auto x0 = node_0[0];
        const auto y0 = node_0[1];
        const auto z0 = node_0[2];
        const auto x1 = node_1[0];
        const auto y1 = node_1[1];
        const auto z1 = node_1[2];
        const auto d = (node_1 - node_0).norm();
        FLOATING_POINT_TYPE analytic_solution = d * (x0 + x1 + 3*y0/2. + 3*y1/2. + 2*z0 + 2*z1 + 5);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);
    }
}
