#pragma once

#include <Eigen/Dense>
#include <Caribou/macros.h>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
#include <Caribou/Geometry/Tetrahedron.h>

TEST(Hexahedron, Linear) {
    using namespace caribou;
    using Hexahedron = caribou::geometry::Hexahedron<Linear>;
    using Tetrahedron = caribou::geometry::Tetrahedron<Linear>;
    using WorldCoordinates = Hexahedron::WorldCoordinates;

    {
        Hexahedron h;
        for (unsigned int i = 0; i < h.number_of_nodes(); ++i) {
            const auto n = h.node(i);
            const auto L = h.L(n);
            std::string shapes;
            bool shapes_are_ok = true;
            for (unsigned int j = 0; j < h.number_of_nodes(); ++j) {
                shapes += std::to_string(L[j]) + " ";
                if (i == j and not IN_CLOSED_INTERVAL(1-EPSILON, L[j], 1+EPSILON)) {
                    shapes_are_ok = false;
                }
                if (i != j and not IN_CLOSED_INTERVAL(0-EPSILON, L[j], 0+EPSILON)) {
                    shapes_are_ok = false;
                }
            }
            EXPECT_TRUE(shapes_are_ok) << "Shape values at node " << i << " are not good: " << shapes;
        }
    }

    {
        Hexahedron h (
            WorldCoordinates(  -1,  1,   -1), WorldCoordinates(   1,  1,   -1),
            WorldCoordinates(   1,  1,    1), WorldCoordinates(  -1,  1,    1),
            WorldCoordinates(-0.5, -1, -0.5), WorldCoordinates( 0.5, -1, -0.5),
            WorldCoordinates( 0.5, -1,  0.5), WorldCoordinates(-0.5, -1,  0.5)
        );

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 8, 1> values;
        values << p1(h.node(0)), p1(h.node(1)), p1(h.node(2)), p1(h.node(3)),
                  p1(h.node(4)), p1(h.node(5)), p1(h.node(6)), p1(h.node(7));
        for (const auto & gauss_node : h.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_FLOAT_EQ(h.interpolate(x, values), p1(h.world_coordinates(x)));
        }

        // Integration
        // Hexa
        FLOATING_POINT_TYPE numerical_solution_hexa = 0;
        for (const auto & gauss_node : h.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = h.jacobian(x).determinant();
            numerical_solution_hexa += p1(h.world_coordinates(x)) * w * detJ;
        }
        // Six tetrahedrons
        FLOATING_POINT_TYPE numerical_solution_tetras = 0;
        Tetrahedron t1 (h.node(0), h.node(5), h.node(1), h.node(6));
        Tetrahedron t2 (h.node(0), h.node(1), h.node(3), h.node(6));
        Tetrahedron t3 (h.node(1), h.node(3), h.node(6), h.node(2));
        Tetrahedron t4 (h.node(6), h.node(3), h.node(0), h.node(7));
        Tetrahedron t5 (h.node(6), h.node(7), h.node(0), h.node(5));
        Tetrahedron t6 (h.node(7), h.node(5), h.node(4), h.node(0));
        for (const auto & t : {t1, t2, t3, t4, t5, t6}) {
            for (const auto & gauss_node : t.gauss_nodes()) {
                const auto x = gauss_node.position;
                const auto w = gauss_node.weight;
                const auto detJ = t.jacobian(x).determinant();
                numerical_solution_tetras += p1(t.world_coordinates(x)) * w * detJ;
            }
        }

        EXPECT_FLOAT_EQ(numerical_solution_hexa, numerical_solution_tetras);
    }
}

TEST(Hexahedron, Quadratic) {
    using namespace caribou;
    using Hexahedron = caribou::geometry::Hexahedron<Quadratic>;
    using Tetrahedron = caribou::geometry::Tetrahedron<Quadratic>;
    {
        Hexahedron h;
        for (unsigned int i = 0; i < h.number_of_nodes(); ++i) {
            const auto n = h.node(i);
            const auto L = h.L(n);
            std::string shapes;
            bool shapes_are_ok = true;
            for (unsigned int j = 0; j < h.number_of_nodes(); ++j) {
                shapes += std::to_string(L[j]) + " ";
                if (i == j and not IN_CLOSED_INTERVAL(1-EPSILON, L[j], 1+EPSILON)) {
                    shapes_are_ok = false;
                }
                if (i != j and not IN_CLOSED_INTERVAL(0-EPSILON, L[j], 0+EPSILON)) {
                    shapes_are_ok = false;
                }
            }
            EXPECT_TRUE(shapes_are_ok) << "Shape values at node " << i << " are not good: " << shapes;
        }
    }

    {
        Hexahedron h;

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 20, 1> values;
        values << p2(h.node( 0)), p2(h.node( 1)), p2(h.node( 2)), p2(h.node( 3)), p2(h.node( 4)), p2(h.node( 5)), p2(h.node( 6)), p2(h.node( 7)), p2(h.node( 8)), p2(h.node( 9)),
                  p2(h.node(10)), p2(h.node(11)), p2(h.node(12)), p2(h.node(13)), p2(h.node(14)), p2(h.node(15)), p2(h.node(16)), p2(h.node(17)), p2(h.node(18)), p2(h.node(19));
        for (const auto & gauss_node : h.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_FLOAT_EQ(h.interpolate(x, values), p2(h.world_coordinates(x)));
        }

        // Integration
        // Hexa
        FLOATING_POINT_TYPE numerical_solution_hexa = 0;
        for (const auto & gauss_node : h.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = h.jacobian(x).determinant();
            numerical_solution_hexa += p2(h.world_coordinates(x)) * w * detJ;
        }
        // Six tetrahedrons
        FLOATING_POINT_TYPE numerical_solution_tetras = 0;
        Tetrahedron t1 ({h.node(0), h.node(5), h.node(1), h.node(6)});
        Tetrahedron t2 ({h.node(0), h.node(1), h.node(3), h.node(6)});
        Tetrahedron t3 ({h.node(1), h.node(3), h.node(6), h.node(2)});
        Tetrahedron t4 ({h.node(6), h.node(3), h.node(0), h.node(7)});
        Tetrahedron t5 ({h.node(6), h.node(7), h.node(0), h.node(5)});
        Tetrahedron t6 ({h.node(7), h.node(5), h.node(4), h.node(0)});
        for (const auto & t : {t1, t2, t3, t4, t5, t6}) {
            for (const auto & gauss_node : t.gauss_nodes()) {
                const auto x = gauss_node.position;
                const auto w = gauss_node.weight;
                const auto detJ = t.jacobian(x).determinant();
                numerical_solution_tetras += p2(t.world_coordinates(x)) * w * detJ;
            }
        }

        EXPECT_FLOAT_EQ(numerical_solution_hexa, numerical_solution_tetras);
    }
}

TEST(Hexahedron, RectangularLinear) {
    using namespace caribou;
    using Hexahedron = caribou::geometry::Hexahedron<Linear>;
    using RectangularHexahedron = caribou::geometry::RectangularHexahedron<Linear>;
    using WorldCoordinates = Hexahedron::WorldCoordinates;
    using Rotation = RectangularHexahedron::Rotation;
    using Size = RectangularHexahedron::Size;

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
        using Segment = caribou::geometry::Segment<_3D, Linear>;
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
        using Triangle = caribou::geometry::Triangle<_3D, Linear>;
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

TEST(Hexahedron, RectangularQuadratic) {
    using namespace caribou;
    using Hexahedron = caribou::geometry::Hexahedron<Quadratic>;
    using RectangularHexahedron = caribou::geometry::RectangularHexahedron<Quadratic>;
    using WorldCoordinates = Hexahedron::WorldCoordinates;
    using Rotation = RectangularHexahedron::Rotation;
    using Size = RectangularHexahedron::Size;

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