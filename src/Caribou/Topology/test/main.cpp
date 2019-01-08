#include <gtest/gtest.h>
#include <Caribou/Topology/Engine/Grid/Grid.h>
#include <Caribou/Topology/Engine/Grid/Cell.h>
#include <Caribou/Algebra/Vector.h>

TEST(Topology, Grid1D) {
    using namespace caribou::topology::engine;
    using Grid = Grid<1>;

    Grid grid(0.25 /* anchor_position */, 2 /* subdivision */, 100 /* size */);

    // Cell size test
    ASSERT_FLOAT_EQ(100/2., grid.H());

    // Cell numbering test
    ASSERT_EQ((Grid::CellIndex) 5, grid.cell_index_at((Grid::GridCoordinates) 5));
    ASSERT_EQ((Grid::GridCoordinates) 5, grid.grid_coordinates_at((Grid::CellIndex) 5));

    // Positioning
    ASSERT_FALSE(grid.contains_position({0.24}));
    ASSERT_TRUE(grid.contains_position({0.25}));
    ASSERT_TRUE(grid.contains_position({50}));
    ASSERT_TRUE(grid.contains_position({100.25}));
    ASSERT_FALSE(grid.contains_position({100.26}));

    ASSERT_EQ((Grid::GridCoordinates) 0, grid.grid_coordinates_at(Grid::WorldCoordinates (0.25)));
    ASSERT_EQ((Grid::GridCoordinates) 0, grid.grid_coordinates_at(Grid::WorldCoordinates (50.24)));
    ASSERT_EQ((Grid::GridCoordinates) 1, grid.grid_coordinates_at(Grid::WorldCoordinates (50.25)));
    ASSERT_EQ((Grid::GridCoordinates) 1, grid.grid_coordinates_at(Grid::WorldCoordinates (100)));
}

TEST(Topology, Grid2D) {
    using namespace caribou::topology::engine;
    using Grid = Grid<2>;

    Grid grid({0.25, 0.5} /* anchor_position */, {4, 4} /* subdivision */, {100, 100} /* size */);

    // Cell numbering test
    ASSERT_EQ((Grid::CellIndex) 5, grid.cell_index_at(Grid::GridCoordinates (1,1)));
    ASSERT_EQ(Grid::GridCoordinates (1,1), grid.grid_coordinates_at((Grid::CellIndex) 5));

    // Positioning
    ASSERT_FALSE(grid.contains_position({0.24,0.5}));
    ASSERT_FALSE(grid.contains_position({0.25,0.49}));
    ASSERT_TRUE(grid.contains_position({0.25,0.5}));
    ASSERT_TRUE(grid.contains_position({50,50}));
    ASSERT_TRUE(grid.contains_position({100.25,100.5}));
    ASSERT_FALSE(grid.contains_position({100.25,100.51}));
    ASSERT_FALSE(grid.contains_position({100.26,100.5}));

}

TEST(Topology, Grid3D) {
    using namespace caribou::topology::engine;
    using Grid = Grid<3>;

    Grid grid({0.25, 0.5, 0.75} /* anchor_position */, {2, 2, 2} /* subdivision */, {100, 100, 100} /* size */);

    // Cell numbering test
    ASSERT_EQ((Grid::CellIndex) 5, grid.cell_index_at(Grid::GridCoordinates (1,0,1)));
    ASSERT_EQ(Grid::GridCoordinates (1,0,1), grid.grid_coordinates_at((Grid::CellIndex) 5));

    // Positioning
    ASSERT_FALSE(grid.contains_position({0.24,0.50, 0.75}));
    ASSERT_FALSE(grid.contains_position({0.25,0.49,0.75}));
    ASSERT_FALSE(grid.contains_position({0.25,0.50,0.74}));
    ASSERT_TRUE(grid.contains_position({0.25,0.50,0.75}));
    ASSERT_TRUE(grid.contains_position({50,50,50}));
    ASSERT_TRUE(grid.contains_position({100.25,100.50,100.75}));
    ASSERT_FALSE(grid.contains_position({100.26,100.50,100.75}));
    ASSERT_FALSE(grid.contains_position({100.25,100.51,100.75}));
    ASSERT_FALSE(grid.contains_position({100.25,100.50,100.76}));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
