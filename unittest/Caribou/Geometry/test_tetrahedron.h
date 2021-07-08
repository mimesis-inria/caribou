#pragma once

#include <Eigen/Dense>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Tetrahedron.h>

TEST(Tetrahedron, Linear) {
    using namespace caribou;

    {
        using Tetrahedron = caribou::geometry::Tetrahedron<Linear>;
        using Triangle = caribou::geometry::Triangle<3, Linear>;
        using Edge = caribou::geometry::Segment<3, Linear>;
        using LocalCoordinates = Tetrahedron::LocalCoordinates;
        using WorldCoordinates = Tetrahedron::WorldCoordinates;

        Tetrahedron t (
            WorldCoordinates({50, 50, 0}),
            WorldCoordinates({60, 50, 0}),
            WorldCoordinates({55, 55, 0}),
            WorldCoordinates({55, 52.5, -5})
        );

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 4, 1> values (p1(t.node(0)), p1(t.node(1)), p1(t.node(2)), p1(t.node(3)));
        for (const auto & gauss_node : t.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(t.interpolate(x, values), p1(t.world_coordinates(x)));
        }

        // Inverse transformation
        for (const auto & gauss_node : t.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, t.local_coordinates(t.world_coordinates(gauss_node.position)), 1e-5);
        }

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                t.local_coordinates(t.center())
            };

            for (UNSIGNED_INTEGER_TYPE face_id = 0; face_id < t.number_of_boundary_elements(); ++face_id) {
                Triangle face = t.boundary_element(face_id);
                for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < face.number_of_gauss_nodes();++gauss_node_id) {
                    inside_points.push_back(t.local_coordinates(face.world_coordinates(face.gauss_node(gauss_node_id).position)));
                }
                for (UNSIGNED_INTEGER_TYPE edge_id = 0; edge_id < face.number_of_boundary_elements(); ++edge_id) {
                    Edge edge = face.boundary_element(edge_id);
                    for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < edge.number_of_gauss_nodes();++gauss_node_id) {
                        inside_points.push_back(t.local_coordinates(face.world_coordinates(face.gauss_node(gauss_node_id).position)));
                    }
                }
            }

            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < t.number_of_nodes();++node_id) {
                inside_points.push_back(t.local_coordinates(t.node(node_id)));
            }
            for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < t.number_of_gauss_nodes();++gauss_node_id) {
                inside_points.push_back(t.gauss_node(gauss_node_id).position);
            }
            for (const auto & p : inside_points) {
                ASSERT_TRUE(t.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << ", " << p[2] << "] is found outside the element, but it should be inside.";
            }

            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {};
            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < t.number_of_nodes();++node_id) {
                const WorldCoordinates center_to_point = t.node(node_id) - t.center();
                const auto length = center_to_point.norm();
                const auto n = center_to_point / length;
                outside_points.emplace_back(t.center() + n*length*1.01);
            }
            for (UNSIGNED_INTEGER_TYPE face_id = 0; face_id < t.number_of_boundary_elements(); ++face_id) {
                Triangle face = t.boundary_element(face_id);
                for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < face.number_of_gauss_nodes(); ++gauss_node_id) {
                    const WorldCoordinates p = face.world_coordinates(face.gauss_node(gauss_node_id).position);
                    const WorldCoordinates center_to_point = p - t.center();
                    outside_points.push_back(t.local_coordinates(t.center() + center_to_point*1.01));
                }
                for (UNSIGNED_INTEGER_TYPE edge_id = 0; edge_id < face.number_of_boundary_elements(); ++edge_id) {
                    Edge edge = face.boundary_element(edge_id);
                    for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < edge.number_of_gauss_nodes();++gauss_node_id) {
                        const WorldCoordinates p = edge.world_coordinates(edge.gauss_node(gauss_node_id).position);
                        const WorldCoordinates center_to_point = p - t.center();
                        outside_points.push_back(t.local_coordinates(t.center() + center_to_point*1.01));
                    }
                }
            }
            for (const auto & p : outside_points) {
                ASSERT_FALSE(t.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << ", " << p[2] << "] is found inside the element, but it should be outside.";
            }
        }
    }
}

TEST(Tetrahedron, Quadratic) {
    using namespace caribou;

    {
        using Tetrahedron = caribou::geometry::Tetrahedron<Quadratic>;
        using Triangle = caribou::geometry::Triangle<3, Quadratic>;
        using Edge = caribou::geometry::Segment<3, Quadratic>;
        using LocalCoordinates = Tetrahedron::LocalCoordinates;
        using WorldCoordinates = Tetrahedron::WorldCoordinates;

        Tetrahedron t (
            WorldCoordinates({50, 50, 0}),
            WorldCoordinates({60, 50, 0}),
            WorldCoordinates({55, 55, 0}),
            WorldCoordinates({55, 52.5, -5})
        );

        // Interpolation
        Eigen::Matrix<FLOATING_POINT_TYPE, 10, 1> values;
        values << p2(t.node(0)), p2(t.node(1)), p2(t.node(2)), p2(t.node(3)), p2(t.node(4)),
                  p2(t.node(5)), p2(t.node(6)), p2(t.node(7)), p2(t.node(8)), p2(t.node(9));
        for (const auto & gauss_node : t.gauss_nodes()) {
            const auto x = gauss_node.position;
            EXPECT_DOUBLE_EQ(t.interpolate(x, values), p2(t.world_coordinates(x)));
        }

        // Inverse transformation
        for (const auto & gauss_node : t.gauss_nodes()) {
            EXPECT_MATRIX_NEAR(gauss_node.position, t.local_coordinates(t.world_coordinates(gauss_node.position)), 1e-5);
        }

        // Contains point
        {
            FLOATING_POINT_TYPE epsilon = 0.000001;
            // Test that all the following points are INSIDE the element
            std::vector<LocalCoordinates> inside_points = {
                t.local_coordinates(t.center())
            };

            for (UNSIGNED_INTEGER_TYPE face_id = 0; face_id < t.number_of_boundary_elements(); ++face_id) {
                Triangle face = t.boundary_element(face_id);
                for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < face.number_of_gauss_nodes(); ++gauss_node_id) {
                    inside_points.push_back(t.local_coordinates(face.world_coordinates(face.gauss_node(gauss_node_id).position)));
                }
                for (UNSIGNED_INTEGER_TYPE edge_id = 0; edge_id < face.number_of_boundary_elements(); ++edge_id) {
                    Edge edge = face.boundary_element(edge_id);
                    for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < edge.number_of_gauss_nodes();++gauss_node_id) {
                        inside_points.push_back(t.local_coordinates(face.world_coordinates(face.gauss_node(gauss_node_id).position)));
                    }
                }
            }

            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < t.number_of_nodes();++node_id) {
                inside_points.push_back(t.local_coordinates(t.node(node_id)));
            }
            for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < t.number_of_gauss_nodes();++gauss_node_id) {
                inside_points.push_back(t.gauss_node(gauss_node_id).position);
            }
            for (const auto & p : inside_points) {
                ASSERT_TRUE(t.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << ", " << p[2] << "] is found outside the element, but it should be inside.";
            }

            // Test that all the following points are OUTSIDE the element
            std::vector<LocalCoordinates> outside_points = {};
            for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < t.number_of_nodes();++node_id) {
                const WorldCoordinates center_to_point = t.node(node_id) - t.center();
                const auto length = center_to_point.norm();
                const auto n = center_to_point / length;
                outside_points.emplace_back(t.center() + n*length*1.01);
            }
            for (UNSIGNED_INTEGER_TYPE face_id = 0; face_id < t.number_of_boundary_elements(); ++face_id) {
                Triangle face = t.boundary_element(face_id);
                for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < face.number_of_gauss_nodes(); ++gauss_node_id) {
                    const WorldCoordinates p = face.world_coordinates(face.gauss_node(gauss_node_id).position);
                    const WorldCoordinates center_to_point = p - t.center();
                    outside_points.push_back(t.local_coordinates(t.center() + center_to_point*1.01));
                }
                for (UNSIGNED_INTEGER_TYPE edge_id = 0; edge_id < face.number_of_boundary_elements(); ++edge_id) {
                    Edge edge = face.boundary_element(edge_id);
                    for (UNSIGNED_INTEGER_TYPE gauss_node_id = 0; gauss_node_id < edge.number_of_gauss_nodes();++gauss_node_id) {
                        const WorldCoordinates p = edge.world_coordinates(edge.gauss_node(gauss_node_id).position);
                        const WorldCoordinates center_to_point = p - t.center();
                        outside_points.push_back(t.local_coordinates(t.center() + center_to_point*1.01));
                    }
                }
            }
            for (const auto & p : outside_points) {
                ASSERT_FALSE(t.contains_local(p, epsilon)) <<
                "Local point [" << p[0] << ", " << p[1] << ", " << p[2] << "] is found inside the element, but it should be outside.";
            }
        }
    }
}