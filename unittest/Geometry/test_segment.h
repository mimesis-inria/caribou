#pragma once

#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Segment.h>

TEST(Segment, Linear) {
    using namespace caribou;

    // Shape functions
    {
        using Segment = caribou::geometry::Segment<_3D, Linear>;
        using LocalCoordinates = Segment::LocalCoordinates ;
        Segment s;
        EXPECT_EQ(s.L(LocalCoordinates(-1))[0], 1);
        EXPECT_EQ(s.L(LocalCoordinates(-1))[1], 0);
        EXPECT_EQ(s.L(LocalCoordinates(+1))[0], 0);
        EXPECT_EQ(s.L(LocalCoordinates(+1))[1], 1);
    }

// 1D
    {
        using Segment = caribou::geometry::Segment<_1D, Linear>;
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
        using Segment = caribou::geometry::Segment<_2D, Linear>;
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
        using Segment = caribou::geometry::Segment<_3D, Linear>;
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

TEST(Segment, Quadratic) {
    using namespace caribou;

    // Shape functions
    {
        using Segment = caribou::geometry::Segment<_3D, Quadratic>;
        using LocalCoordinates = Segment::LocalCoordinates ;
        Segment s;
        EXPECT_EQ(s.L(LocalCoordinates(-1))[0], 1);
        EXPECT_EQ(s.L(LocalCoordinates(-1))[1], 0);
        EXPECT_EQ(s.L(LocalCoordinates(-1))[2], 0);
        EXPECT_EQ(s.L(LocalCoordinates(+1))[0], 0);
        EXPECT_EQ(s.L(LocalCoordinates(+1))[1], 1);
        EXPECT_EQ(s.L(LocalCoordinates(+1))[2], 0);
        EXPECT_EQ(s.L(LocalCoordinates( 0))[0], 0);
        EXPECT_EQ(s.L(LocalCoordinates( 0))[1], 0);
        EXPECT_EQ(s.L(LocalCoordinates( 0))[2], 1);
    }

// 1D
    {
        using Segment = caribou::geometry::Segment<_1D, Quadratic>;
        using WorldCoordinates = Segment::WorldCoordinates;
        Segment segment(WorldCoordinates(-5.5), WorldCoordinates(1.1), WorldCoordinates(-2.2));

        // Center
        EXPECT_DOUBLE_EQ(segment.center()[0], -2.2);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p2(segment.node(0)), p2(segment.node(1)), p2(segment.node(2)));
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(segment.interpolate(x, values), p2(segment.world_coordinates(x)));
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = segment.jacobian(x)[0];
            numerical_solution += p2(segment.world_coordinates(x)) * w * detJ;
        }
        FLOATING_POINT_TYPE x0 = segment.node(0)[0];
        FLOATING_POINT_TYPE x1 = segment.node(1)[0];
        auto analytic_solution = static_cast<FLOATING_POINT_TYPE>((2 * x1*x1*x1 / 3. + 5*x1) - (2 * x0*x0*x0 / 3. + 5*x0)) ;
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);
    }

// 2D
    {
        using Segment = caribou::geometry::Segment<_2D, Quadratic>;
        using WorldCoordinates  = Segment::WorldCoordinates ;

        WorldCoordinates node_0 {-1.5, -1.5};
        WorldCoordinates node_1 {5.5, 5.5};
        WorldCoordinates node_2 {2, 2};

        Segment segment(node_0, node_1, node_2);

        // Center
        WorldCoordinates center_node = node_0 + (node_1 - node_0).normalized()*(node_1-node_0).norm()/2.;
        EXPECT_DOUBLE_EQ(segment.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(segment.center()[1], center_node[1]);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p2(node_0), p2(node_1), p2(node_2));
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(segment.interpolate(x, values), p2(segment.world_coordinates(x)));
        }

        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  segment.jacobian(x);
            const auto detJ = J.norm();
            numerical_solution += p2(segment.world_coordinates(x)) * w * detJ;
        }

        // Integration
        const auto x0 = node_0[0]; // a
        const auto y0 = node_0[1]; // c
        const auto x1 = node_1[0]; // b
        const auto y1 = node_1[1]; // d
        const auto d = (node_1 - node_0).norm();
        FLOATING_POINT_TYPE analytic_solution = d * (1/3.*x0*(2*y0+y1) + 1/3.*x1*(2*y1+y0) + y0*y0 + y0*y1 + y1*y1 + 5);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);
    }

// 3D
    {
        using Segment = caribou::geometry::Segment<_3D, Quadratic>;
        using WorldCoordinates  = Segment::WorldCoordinates ;

        WorldCoordinates node_0 {-1.5, -1.5, -5.2};
        WorldCoordinates node_1 {5.5, 5.5, 54.3};
        WorldCoordinates node_2 {2, 2, 24.55};

        // Center
        Segment segment(node_0, node_1, node_2);
        WorldCoordinates center_node = node_0 + (node_1 - node_0).normalized()*(node_1-node_0).norm()/2.;
        EXPECT_DOUBLE_EQ(segment.center()[0], center_node[0]);
        EXPECT_DOUBLE_EQ(segment.center()[1], center_node[1]);
        EXPECT_DOUBLE_EQ(segment.center()[2], center_node[2]);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1> values (p2(node_0), p2(node_1), p2(node_2));
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(segment.interpolate(x, values), p2(segment.world_coordinates(x)));
        }

        // Integration
        FLOATING_POINT_TYPE numerical_solution = 0;
        for (const auto & gauss_node : segment.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J =  segment.jacobian(x);
            const auto detJ = J.norm();
            numerical_solution += p2(segment.world_coordinates(x)) * w * detJ;
        }

        const auto a = node_0[0]; // x0
        const auto b = node_1[0]; // x1
        const auto c = node_0[1]; // y0
        const auto d = node_1[1]; // y1
        const auto e = node_0[2]; // z0
        const auto f = node_1[2]; // z1
        const auto n = (node_1 - node_0).norm();
        FLOATING_POINT_TYPE analytic_solution = n * ((30 + 2*a*(2*c + d) + 2*b*(c + 2*d) + 6*c*e + 3*d*e + 8*e*e + 3*c*f + 6*d*f + 8*e*f + 8*f*f)/6);
        EXPECT_DOUBLE_EQ(numerical_solution, analytic_solution);
    }
}
