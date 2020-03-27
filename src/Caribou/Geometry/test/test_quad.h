#pragma once

#include <Eigen/Dense>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Quad.h>

TEST(Quad, Linear) {
    using namespace caribou;

    // 2D
    {
        using Triangle = caribou::geometry::Triangle<_2D, Linear>;
        using Quad = caribou::geometry::Quad<_2D, Linear>;

        using WorldCoordinates = Quad::WorldCoordinates;

        Quad q (
            WorldCoordinates(-5, -5), WorldCoordinates(+5, -5),
            WorldCoordinates(+10, +5), WorldCoordinates(-10, +5)
        );

        EXPECT_EQ(q.number_of_boundary_elements(), 4);

        // Center
        EXPECT_FLOAT_EQ(q.center()[0], 0);
        EXPECT_FLOAT_EQ(q.center()[1], 0);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 4, 1> values (p1(q.node(0)), p1(q.node(1)), p1(q.node(2)), p1(q.node(3)));
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_FLOAT_EQ(q.interpolate(x, values), p1(q.T(x)));
        }

        // Integration
        // Quad
        FLOATING_POINT_TYPE numerical_solution_quad = 0;
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = q.jacobian(x).determinant();
            numerical_solution_quad += p1(q.T(x)) * w * detJ;
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
                numerical_solution_triangles += p1(t.T(x)) * w * detJ;
            }
        }

        EXPECT_FLOAT_EQ(numerical_solution_quad, numerical_solution_triangles);
    }

    // 3D
    {
        using Triangle = caribou::geometry::Triangle<_3D, Linear>;
        using Quad = caribou::geometry::Quad<_3D, Linear>;

        using WorldCoordinates = Quad::WorldCoordinates;

        Quad q (
            WorldCoordinates(-5, -53./15, -53./15), WorldCoordinates(+5, -53./15, -53./15),
            WorldCoordinates(+10,+53./15, +53./15), WorldCoordinates(-10,+53./15, +53./15)
        );

        EXPECT_EQ(q.number_of_boundary_elements(), 4);

        // Center
        EXPECT_FLOAT_EQ(q.center()[0], 0);
        EXPECT_FLOAT_EQ(q.center()[1], 0);
        EXPECT_FLOAT_EQ(q.center()[2], 0);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 4, 1> values (p1(q.node(0)), p1(q.node(1)), p1(q.node(2)), p1(q.node(3)));
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_FLOAT_EQ(q.interpolate(x, values), p1(q.T(x)));
        }

        // Integration
        // Quad
        FLOATING_POINT_TYPE numerical_solution_quad = 0;
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J = q.jacobian(x);
            const auto detJ = J.col(0).cross(J.col(1)).norm();
            numerical_solution_quad += p1(q.T(x)) * w * detJ;
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
                numerical_solution_triangles += p1(t.T(x)) * w * detJ;
            }
        }

        EXPECT_FLOAT_EQ(numerical_solution_quad, numerical_solution_triangles);
    }
}

TEST(Quad, Quadratic) {
    using namespace caribou;

    // 2D
    {
        using Triangle = caribou::geometry::Triangle<_2D, Quadratic>;
        using Quad = caribou::geometry::Quad<_2D, Quadratic>;

        using WorldCoordinates = Quad::WorldCoordinates;

        Quad q ({
            WorldCoordinates(-5, -5), WorldCoordinates(+5, -5),
            WorldCoordinates(+10, +5), WorldCoordinates(-10, +5)
        });

        EXPECT_EQ(q.number_of_boundary_elements(), 4);

        // Center
        EXPECT_FLOAT_EQ(q.center()[0], 0);
        EXPECT_FLOAT_EQ(q.center()[1], 0);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 8, 1> values;
        values << p2(q.node(0)), p2(q.node(1)), p2(q.node(2)), p2(q.node(3)),
                  p2(q.node(4)), p2(q.node(5)), p2(q.node(6)), p2(q.node(7));
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_FLOAT_EQ(q.interpolate(x, values), p2(q.T(x)));
        }

        // Integration
        // Quad
        FLOATING_POINT_TYPE numerical_solution_quad = 0;
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = q.jacobian(x).determinant();
            numerical_solution_quad += p2(q.T(x)) * w * detJ;
        }

        // Two triangles
        FLOATING_POINT_TYPE numerical_solution_triangles = 0;
        Triangle t1 ({q.node(0), q.node(1), q.node(3)});
        Triangle t2 ({q.node(1), q.node(2), q.node(3)});
        for (const auto & t : {t1, t2}) {
            for (const auto & gauss_node : t.gauss_nodes()) {
                const auto x = gauss_node.position;
                const auto w = gauss_node.weight;
                const auto detJ = t.jacobian(x).determinant();
                numerical_solution_triangles += p2(t.T(x)) * w * detJ;
            }
        }

        EXPECT_FLOAT_EQ(numerical_solution_quad, numerical_solution_triangles);
    }

    // 3D
    {
        using Triangle = caribou::geometry::Triangle<_3D, Quadratic>;
        using Quad = caribou::geometry::Quad<_3D, Quadratic>;

        using WorldCoordinates = Quad::WorldCoordinates;

        Quad q ({
            WorldCoordinates(-5, -53./15, -53./15), WorldCoordinates(+5, -53./15, -53./15),
            WorldCoordinates(+10,+53./15, +53./15), WorldCoordinates(-10,+53./15, +53./15)
        });

        EXPECT_EQ(q.number_of_boundary_elements(), 4);

        // Center
        EXPECT_FLOAT_EQ(q.center()[0], 0);
        EXPECT_FLOAT_EQ(q.center()[1], 0);
        EXPECT_FLOAT_EQ(q.center()[2], 0);

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 8, 1> values;
        values << p2(q.node(0)), p2(q.node(1)), p2(q.node(2)), p2(q.node(3)),
                  p2(q.node(4)), p2(q.node(5)), p2(q.node(6)), p2(q.node(7));
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_FLOAT_EQ(q.interpolate(x, values), p2(q.T(x)));
        }

        // Integration
        // Quad
        FLOATING_POINT_TYPE numerical_solution_quad = 0;
        for (const auto & gauss_node : q.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto J = q.jacobian(x);
            const auto detJ = J.col(0).cross(J.col(1)).norm();
            numerical_solution_quad += p2(q.T(x)) * w * detJ;
        }

        // Two triangles
        FLOATING_POINT_TYPE numerical_solution_triangles = 0;
        Triangle t1 ({q.node(0), q.node(1), q.node(3)});
        Triangle t2 ({q.node(1), q.node(2), q.node(3)});
        for (const auto & t : {t1, t2}) {
            for (const auto & gauss_node : t.gauss_nodes()) {
                const auto x = gauss_node.position;
                const auto w = gauss_node.weight;
                const auto J = t.jacobian(x);
                const auto detJ = J.col(0).cross(J.col(1)).norm();
                numerical_solution_triangles += p2(t.T(x)) * w * detJ;
            }
        }

        EXPECT_FLOAT_EQ(numerical_solution_quad, numerical_solution_triangles);
    }
}