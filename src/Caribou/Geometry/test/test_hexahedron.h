#pragma once

#include <Eigen/Dense>
#include <Caribou/macros.h>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Hexahedron.h>
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
            EXPECT_FLOAT_EQ(h.interpolate(x, values), p1(h.T(x)));
        }

        // Integration
        // Hexa
        FLOATING_POINT_TYPE numerical_solution_hexa = 0;
        for (const auto & gauss_node : h.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = h.jacobian(x).determinant();
            numerical_solution_hexa += p1(h.T(x)) * w * detJ;
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
                numerical_solution_tetras += p1(t.T(x)) * w * detJ;
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
            EXPECT_FLOAT_EQ(h.interpolate(x, values), p2(h.T(x)));
        }

        // Integration
        // Hexa
        FLOATING_POINT_TYPE numerical_solution_hexa = 0;
        for (const auto & gauss_node : h.gauss_nodes()) {
            const auto x = gauss_node.position;
            const auto w = gauss_node.weight;
            const auto detJ = h.jacobian(x).determinant();
            numerical_solution_hexa += p2(h.T(x)) * w * detJ;
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
                numerical_solution_tetras += p2(t.T(x)) * w * detJ;
            }
        }

        EXPECT_FLOAT_EQ(numerical_solution_hexa, numerical_solution_tetras);
    }
}