#include <gtest/gtest.h>
#include <Caribou/Topology/Engine/Grid/Grid.h>
#include <Caribou/Topology/Engine/Grid/GridContainer.h>
#include <Caribou/Topology/Engine/Grid/Cell.h>
#include <Caribou/Algebra/Vector.h>

template<typename T>
bool list_are_equals (const std::list<T> & l1, const std::list<T> & l2)
{
    if (l1.size() != l2.size())
        return false;
    for (auto it1 = l1.begin(), it2 = l2.begin(); it1 != l1.end(); ++it1, ++it2)
        if (*it1 != *it2)
            return false;
    return true;
}

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

    // Cells queries
    ASSERT_TRUE(list_are_equals({}, grid.cells_enclosing(-50, 0, 0.24)));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_enclosing(-50, 25, 50.24)));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_enclosing(0.25, 25, 50.24)));
    ASSERT_TRUE(list_are_equals({1}, grid.cells_enclosing(50.25, 75, 100.24)));
    ASSERT_TRUE(list_are_equals({0,1}, grid.cells_enclosing(25, 75)));
    ASSERT_TRUE(list_are_equals({0,1}, grid.cells_enclosing(0.25, 101)));
}

TEST(Topology, GridContainer1D) {
    using namespace caribou::topology::engine;
    using Cell = Cell<1>;
    using Grid = GridContainer<1, Cell>;

    Grid grid(0.25 /* anchor_position */, 2 /* subdivision */, 100 /* size */);
}

TEST(Topology, Grid2D) {
    using namespace caribou::topology::engine;
    using Grid = Grid<2>;

    Grid grid({0.25, 0.5} /* anchor_position */, {2, 2} /* subdivision */, {100, 100} /* size */);

    // Cell numbering test
    ASSERT_EQ((Grid::CellIndex) 3, grid.cell_index_at(Grid::GridCoordinates (1,1)));
    ASSERT_EQ(Grid::GridCoordinates (1,1), grid.grid_coordinates_at((Grid::CellIndex) 3));

    // Positioning
    ASSERT_FALSE(grid.contains_position({0.24,0.5}));
    ASSERT_FALSE(grid.contains_position({0.25,0.49}));
    ASSERT_TRUE(grid.contains_position({0.25,0.5}));
    ASSERT_TRUE(grid.contains_position({50,50}));
    ASSERT_TRUE(grid.contains_position({100.25,100.5}));
    ASSERT_FALSE(grid.contains_position({100.25,100.51}));
    ASSERT_FALSE(grid.contains_position({100.26,100.5}));

    // Cells queries
    ASSERT_TRUE(list_are_equals({0,2}, grid.cells_enclosing(
            Grid::WorldCoordinates ({-50, 50}),
            Grid::WorldCoordinates ({25, 25}),
            Grid::WorldCoordinates ({25, 75})
    )));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_enclosing(
            Grid::WorldCoordinates ({0.25, 0.5}),
            Grid::WorldCoordinates ({25, 25}),
            Grid::WorldCoordinates ({50.24, 50.49})
            )));
    ASSERT_TRUE(list_are_equals({1}, grid.cells_enclosing(
            Grid::WorldCoordinates ({50.25, 0.5}),
            Grid::WorldCoordinates ({75, 25}),
            Grid::WorldCoordinates ({100.24, 50.49})
    )));
    ASSERT_TRUE(list_are_equals({0,1,2,3}, grid.cells_enclosing(
            Grid::WorldCoordinates ({0.25, 0.5}),
            Grid::WorldCoordinates ({75.25, 75.50})
    )));

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

    // Cells queries
    ASSERT_TRUE(list_are_equals({0,2,4,6}, grid.cells_enclosing(
            Grid::WorldCoordinates ({-50, 50, 50}),
            Grid::WorldCoordinates ({25, 25, 25}),
            Grid::WorldCoordinates ({25, 75, 75})
    )));

    ASSERT_TRUE(list_are_equals({0}, grid.cells_enclosing(
            Grid::WorldCoordinates ({0.25, 0.5, 0.75}),
            Grid::WorldCoordinates ({25, 25, 25}),
            Grid::WorldCoordinates ({50.24, 50.49, 50.74})
    )));
    ASSERT_TRUE(list_are_equals({7}, grid.cells_enclosing(
            Grid::WorldCoordinates ({50.25, 50.5, 50.75}),
            Grid::WorldCoordinates ({75, 75, 75}),
            Grid::WorldCoordinates ({100.24, 100.49, 100.74})
    )));
    ASSERT_TRUE(list_are_equals({0,1,2,3,4,5,6,7}, grid.cells_enclosing(
            Grid::WorldCoordinates ({0.25, 0.5, 0.75}),
            Grid::WorldCoordinates ({100.24, 100.49, 100.74})
    )));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
