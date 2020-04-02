#pragma once

#include <Eigen/Dense>
#include <Caribou/Constants.h>
#include <Caribou/Geometry/Tetrahedron.h>

TEST(Tetrahedron, Linear) {
    using namespace caribou;

    {
        using Tetrahedron = caribou::geometry::Tetrahedron<Linear>;

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
            EXPECT_FLOAT_EQ(t.interpolate(x, values), p1(t.world_coordinates(x)));
        }
    }
}

TEST(Tetrahedron, Quadratic) {
    using namespace caribou;

    {
        using Tetrahedron = caribou::geometry::Tetrahedron<Quadratic>;

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
            EXPECT_FLOAT_EQ(t.interpolate(x, values), p2(t.world_coordinates(x)));
        }
    }
}