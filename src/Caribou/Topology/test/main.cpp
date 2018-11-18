#include <gtest/gtest.h>
#include <Caribou/Topology/Engine/Grid/Grid.inl>
#include <Caribou/Topology/Engine/Grid/Cell.h>

TEST(Topology, Grid2D) {
    using namespace caribou::topology::engine;
    using Cell = Cell<2>;
    using Grid = Grid<Cell>;
    using Coord = Grid::VecFloat;

    Grid grid(
            {0.75, 0.5} /* anchor */,
            {4, 4} /* subdivision */,
            {100, 100} /* size */
    );

    // Grid's cell accessor test
    ASSERT_THROW(grid.get({4,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,4}), std::out_of_range);

    // Cell size test
    const Cell & cell = grid.get({0,0});
    ASSERT_FLOAT_EQ(100/4., grid.cell_size(cell)[0]);

    // Cell numbering test
    ASSERT_EQ((Grid::Index)0, grid.get({0,0}).index());
    ASSERT_EQ((Grid::Index)15, grid.get({3,3}).index());

    // Node numbering test
    ASSERT_EQ((Grid::Index)24, grid.nodes({3,3})[2]);
    ASSERT_EQ(Grid::VecInt ({2,1}), grid.grid_coordinates(grid.cell_index({2,1})));

    // Positioning test
    Grid::Index node_id = grid.nodes({1,2})[0];
    ASSERT_EQ((Grid::Index)11, node_id);
    ASSERT_EQ(Grid::VecFloat({25.75, 50.5}), grid.position(node_id));

    // Subdivision test
    Cell & outer_cell = grid.get({2,1});
    outer_cell.subdivide();
    Cell & inner_cell = outer_cell.child({1,0});

    std::array<Coord, 4> p = grid.positions(inner_cell);

    ASSERT_EQ(p[0], Coord({    63.25,      25.5         }));
    ASSERT_EQ(p[1], Coord({    63.25+12.5, 25.5         }));
    ASSERT_EQ(p[2], Coord({    63.25+12.5, 25.5+12.5    }));
    ASSERT_EQ(p[3], Coord({    63.25,      25.5+12.5    }));
}

TEST(Topology, Grid3D) {
    using namespace caribou::topology::engine;
    using Cell = Cell<3>;
    using Grid = Grid<Cell>;
    using Coord = Grid::VecFloat;

    Grid grid(
            {0.25, 0.5, 0.75} /* anchor */,
            {2, 2, 2} /* subdivision */,
            {100, 100, 100} /* size */
    );

    // Grid's cell accessor test
    ASSERT_THROW(grid.get({2,0,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,2,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,0,2}), std::out_of_range);

    // Cell size test
    const Cell & cell = grid.get({0,0,0});
    ASSERT_FLOAT_EQ(100/2., grid.cell_size(cell)[0]);

    // Cell numbering test
    ASSERT_EQ((Grid::Index)0, grid.get({0,0,0}).index());
    ASSERT_EQ((Grid::Index)3, grid.get({1,1,0}).index());
    ASSERT_EQ((Grid::Index)6, grid.get({0,1,1}).index());

    // Node numbering test
    Grid::Index node_id = grid.nodes({1,0,1})[0];
    ASSERT_EQ((Grid::Index)10, node_id);

    // Positioning test
    ASSERT_EQ(Grid::VecFloat({50.25, 0.5, 50.75}), grid.position(node_id));

    // Subdivision test
    Cell & outer_cell = grid.get({1,0,1});
    outer_cell.subdivide();
    Cell & inner_cell = outer_cell.child({0,1,0});

    std::array<Coord, 8> p = grid.positions(inner_cell);

    ASSERT_EQ(p[0], Coord({    50.25,    25.5,    50.75       }));
    ASSERT_EQ(p[1], Coord({    50.25+25, 25.5,    50.75       }));
    ASSERT_EQ(p[2], Coord({    50.25+25, 25.5+25, 50.75       }));
    ASSERT_EQ(p[3], Coord({    50.25,    25.5+25, 50.75       }));
    ASSERT_EQ(p[4], Coord({    50.25,    25.5,    50.75+25    }));
    ASSERT_EQ(p[5], Coord({    50.25+25, 25.5,    50.75+25    }));
    ASSERT_EQ(p[6], Coord({    50.25+25, 25.5+25, 50.75+25    }));
    ASSERT_EQ(p[7], Coord({    50.25,    25.5+25, 50.75+25    }));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
