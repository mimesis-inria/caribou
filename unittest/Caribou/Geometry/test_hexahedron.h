#pragma once

#include <Eigen/Dense>
#include <Caribou/macros.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Tetrahedron.h>

TEST(Hexahedron, Linear) {
    using namespace caribou;
    using Hexahedron = caribou::geometry::Hexahedron;
    using Tetrahedron = caribou::geometry::Tetrahedron;
    using Quad = caribou::geometry::Quad<3>;
    using Triangle = caribou::geometry::Triangle<3>;
    using Edge = caribou::geometry::Segment<3>;
    using LocalCoordinates = Hexahedron::LocalCoordinates;
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
            EXPECT_NEAR(h.interpolate(x, values), p1(h.world_coordinates(x)), 1e-10);
        }

        // Inverse transformation
        for (const auto & gauss_node : h.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, h.local_coordinates(h.world_coordinates(gauss_node.position)), 1e-5);
        }

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                h.local_coordinates(h.center())
            };

            for (UNSIGNED_INTEGER_TYPE face_id = 0; face_id < h.number_of_boundary_elements(); ++face_id) {
                Quad face = h.boundary_element(face_id);
                for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < face.number_of_gauss_nodes();++gauss_node_id) {
                    inside_points.push_back(h.local_coordinates(face.world_coordinates(face.gauss_node(gauss_node_id).position)));
                }
                for (UNSIGNED_INTEGER_TYPE edge_id = 0; edge_id < face.number_of_boundary_elements(); ++edge_id) {
                    Edge edge = face.boundary_element(edge_id);
                    for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < edge.number_of_gauss_nodes();++gauss_node_id) {
                        inside_points.push_back(h.local_coordinates(face.world_coordinates(face.gauss_node(gauss_node_id).position)));
                    }
                }
            }

            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < h.number_of_nodes();++node_id) {
                inside_points.push_back(h.local_coordinates(h.node(node_id)));
            }
            for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < h.number_of_gauss_nodes();++gauss_node_id) {
                inside_points.push_back(h.gauss_node(gauss_node_id).position);
            }
            for (const auto & p : inside_points) {
                ASSERT_TRUE(h.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << ", " << p[2] << "] is found outside the element, but it should be inside.";
            }

            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {};
            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < h.number_of_nodes();++node_id) {
                const WorldCoordinates center_to_point = h.node(node_id) - h.center();
                outside_points.emplace_back(h.center() + center_to_point*1.01);
            }
            for (UNSIGNED_INTEGER_TYPE face_id = 0; face_id < h.number_of_boundary_elements(); ++face_id) {
                Quad face = h.boundary_element(face_id);
                for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < face.number_of_gauss_nodes(); ++gauss_node_id) {
                    const WorldCoordinates p = face.world_coordinates(face.gauss_node(gauss_node_id).position);
                    const WorldCoordinates center_to_point = p - h.center();
                    outside_points.push_back(h.local_coordinates(h.center() + center_to_point*1.01));
                }
                for (UNSIGNED_INTEGER_TYPE edge_id = 0; edge_id < face.number_of_boundary_elements(); ++edge_id) {
                    Edge edge = face.boundary_element(edge_id);
                    for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < edge.number_of_gauss_nodes();++gauss_node_id) {
                        const WorldCoordinates p = edge.world_coordinates(edge.gauss_node(gauss_node_id).position);
                        const WorldCoordinates center_to_point = p - h.center();
                        outside_points.push_back(h.local_coordinates(h.center() + center_to_point*1.01));
                    }
                }
            }
            for (const auto & p : outside_points) {
                ASSERT_FALSE(h.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << ", " << p[2] << "] is found inside the element, but it should be outside.";
            }
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

        EXPECT_DOUBLE_EQ(numerical_solution_hexa, numerical_solution_tetras);

        // Intersection with triangles (test taken from a failing scene)
        {
            Hexahedron h (
                    WorldCoordinates(194.948, -111.776, -1547.85), WorldCoordinates(195.837, -111.776, -1547.85),
                    WorldCoordinates(195.837, -110.461, -1547.85), WorldCoordinates(194.948, -110.461, -1547.85),
                    WorldCoordinates(194.948, -111.776, -1546.86), WorldCoordinates(195.837, -111.776, -1546.86),
                    WorldCoordinates(195.837, -110.461, -1546.86), WorldCoordinates(194.948, -110.461, -1546.86)
            );

            Triangle t (WorldCoordinates {195.402, -112.092, -1543.81}, WorldCoordinates {196.468, -108.742, -1546.42}, WorldCoordinates {193.869, -110.918, -1548.45});

            ASSERT_TRUE(h.intersects(t, 0));
        }
    }
}
