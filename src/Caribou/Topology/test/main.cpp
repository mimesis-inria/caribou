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

    // General properties
    ASSERT_EQ(grid.number_of_nodes(), 3);

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

    // Node queries
    ASSERT_EQ(grid.node(0), 0.25);
    ASSERT_EQ(grid.node(1), 50.25);
    ASSERT_EQ(grid.node(2), 100.25);

    // Edge queries
    ASSERT_EQ(grid.number_of_edges(), 2);
    ASSERT_EQ(grid.edge(0), Grid::Edge (0.25, 50.25));
    ASSERT_EQ(grid.edge(1), Grid::Edge (50.25, 100.25));
}

TEST(Topology, Grid2D) {
    using namespace caribou::topology::engine;
    using Grid = Grid<2>;

    Grid grid({0.25, 0.5} /* anchor_position */, {2, 2} /* subdivision */, {100, 100} /* size */);

    // General properties
    ASSERT_EQ(grid.number_of_nodes(), 9);

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

    // Node queries
    ASSERT_EQ(grid.node(0), Grid::WorldCoordinates ({0.25, 0.5}));
    ASSERT_EQ(grid.node(1), Grid::WorldCoordinates ({50.25, 0.5}));
    ASSERT_EQ(grid.node(2), Grid::WorldCoordinates ({100.25, 0.5}));
    ASSERT_EQ(grid.node(3), Grid::WorldCoordinates ({0.25, 50.5}));
    ASSERT_EQ(grid.node(4), Grid::WorldCoordinates ({50.25, 50.5}));
    ASSERT_EQ(grid.node(5), Grid::WorldCoordinates ({100.25, 50.5}));
    ASSERT_EQ(grid.node(6), Grid::WorldCoordinates ({0.25, 100.5}));
    ASSERT_EQ(grid.node(7), Grid::WorldCoordinates ({50.25, 100.5}));
    ASSERT_EQ(grid.node(8), Grid::WorldCoordinates ({100.25, 100.5}));

    // Edge queries
    ASSERT_EQ(grid.number_of_edges(), 12);
    ASSERT_EQ(grid.edge(0), Grid::Edge (grid.node(0), grid.node(1)));
    ASSERT_EQ(grid.edge(1), Grid::Edge (grid.node(1), grid.node(2)));
    ASSERT_EQ(grid.edge(2), Grid::Edge (grid.node(0), grid.node(3)));
    ASSERT_EQ(grid.edge(3), Grid::Edge (grid.node(1), grid.node(4)));
    ASSERT_EQ(grid.edge(4), Grid::Edge (grid.node(2), grid.node(5)));
    ASSERT_EQ(grid.edge(5), Grid::Edge (grid.node(3), grid.node(4)));
    ASSERT_EQ(grid.edge(6), Grid::Edge (grid.node(4), grid.node(5)));
    ASSERT_EQ(grid.edge(7), Grid::Edge (grid.node(3), grid.node(6)));
    ASSERT_EQ(grid.edge(8), Grid::Edge (grid.node(4), grid.node(7)));
    ASSERT_EQ(grid.edge(9), Grid::Edge (grid.node(5), grid.node(8)));
    ASSERT_EQ(grid.edge(10), Grid::Edge (grid.node(6), grid.node(7)));
    ASSERT_EQ(grid.edge(11), Grid::Edge (grid.node(7), grid.node(8)));

}

TEST(Topology, Grid3D) {
    using namespace caribou::topology::engine;
    using Grid = Grid<3>;

    Grid grid({0.25, 0.5, 0.75} /* anchor_position */, {2, 2, 2} /* subdivision */, {100, 100, 100} /* size */);
    
    // General properties
    ASSERT_EQ(grid.number_of_nodes(), 27);

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

    // Node queries
    ASSERT_EQ(grid.node(0), Grid::WorldCoordinates ({0.25, 0.5, 0.75}));
    ASSERT_EQ(grid.node(2), Grid::WorldCoordinates ({100.25, 0.5, 0.75}));
    ASSERT_EQ(grid.node(6), Grid::WorldCoordinates ({0.25, 100.5, 0.75}));
    ASSERT_EQ(grid.node(8), Grid::WorldCoordinates ({100.25, 100.5, 0.75}));

    ASSERT_EQ(grid.node(18), Grid::WorldCoordinates ({0.25, 0.5, 100.75}));
    ASSERT_EQ(grid.node(20), Grid::WorldCoordinates ({100.25, 0.5, 100.75}));
    ASSERT_EQ(grid.node(24), Grid::WorldCoordinates ({0.25, 100.5, 100.75}));
    ASSERT_EQ(grid.node(26), Grid::WorldCoordinates ({100.25, 100.5, 100.75}));

    // Edge queries
    ASSERT_EQ(grid.number_of_edges(), 3*12+2*9);
    // First slice (2D grid)
    ASSERT_EQ(grid.edge(0), Grid::Edge (grid.node(0), grid.node(1)));
    ASSERT_EQ(grid.edge(1), Grid::Edge (grid.node(1), grid.node(2)));
    ASSERT_EQ(grid.edge(2), Grid::Edge (grid.node(0), grid.node(3)));
    ASSERT_EQ(grid.edge(3), Grid::Edge (grid.node(1), grid.node(4)));
    ASSERT_EQ(grid.edge(4), Grid::Edge (grid.node(2), grid.node(5)));
    ASSERT_EQ(grid.edge(5), Grid::Edge (grid.node(3), grid.node(4)));
    ASSERT_EQ(grid.edge(6), Grid::Edge (grid.node(4), grid.node(5)));
    ASSERT_EQ(grid.edge(7), Grid::Edge (grid.node(3), grid.node(6)));
    ASSERT_EQ(grid.edge(8), Grid::Edge (grid.node(4), grid.node(7)));
    ASSERT_EQ(grid.edge(9), Grid::Edge (grid.node(5), grid.node(8)));
    ASSERT_EQ(grid.edge(10), Grid::Edge (grid.node(6), grid.node(7)));
    ASSERT_EQ(grid.edge(11), Grid::Edge (grid.node(7), grid.node(8)));
    // Between the slice 1 and slice 2
    ASSERT_EQ(grid.edge(12), Grid::Edge (grid.node(0), grid.node(9)));
    ASSERT_EQ(grid.edge(13), Grid::Edge (grid.node(1), grid.node(10)));
    ASSERT_EQ(grid.edge(14), Grid::Edge (grid.node(2), grid.node(11)));
    ASSERT_EQ(grid.edge(15), Grid::Edge (grid.node(3), grid.node(12)));
    ASSERT_EQ(grid.edge(16), Grid::Edge (grid.node(4), grid.node(13)));
    ASSERT_EQ(grid.edge(17), Grid::Edge (grid.node(5), grid.node(14)));
    ASSERT_EQ(grid.edge(18), Grid::Edge (grid.node(6), grid.node(15)));
    ASSERT_EQ(grid.edge(19), Grid::Edge (grid.node(7), grid.node(16)));
    ASSERT_EQ(grid.edge(20), Grid::Edge (grid.node(8), grid.node(17)));
    // second slice (2D grid)
    ASSERT_EQ(grid.edge(21), Grid::Edge (grid.node(9), grid.node(10)));
    ASSERT_EQ(grid.edge(22), Grid::Edge (grid.node(10), grid.node(11)));
    ASSERT_EQ(grid.edge(23), Grid::Edge (grid.node(9), grid.node(12)));
    ASSERT_EQ(grid.edge(24), Grid::Edge (grid.node(10), grid.node(13)));
    ASSERT_EQ(grid.edge(25), Grid::Edge (grid.node(11), grid.node(14)));
    ASSERT_EQ(grid.edge(26), Grid::Edge (grid.node(12), grid.node(13)));
    ASSERT_EQ(grid.edge(27), Grid::Edge (grid.node(13), grid.node(14)));
    ASSERT_EQ(grid.edge(28), Grid::Edge (grid.node(12), grid.node(15)));
    ASSERT_EQ(grid.edge(29), Grid::Edge (grid.node(13), grid.node(16)));
    ASSERT_EQ(grid.edge(30), Grid::Edge (grid.node(14), grid.node(17)));
    ASSERT_EQ(grid.edge(31), Grid::Edge (grid.node(15), grid.node(16)));
    ASSERT_EQ(grid.edge(32), Grid::Edge (grid.node(16), grid.node(17)));
    // Between the slice 2 and slice 3
    ASSERT_EQ(grid.edge(33), Grid::Edge (grid.node(9), grid.node(18)));
    ASSERT_EQ(grid.edge(34), Grid::Edge (grid.node(10), grid.node(19)));
    ASSERT_EQ(grid.edge(35), Grid::Edge (grid.node(11), grid.node(20)));
    ASSERT_EQ(grid.edge(36), Grid::Edge (grid.node(12), grid.node(21)));
    ASSERT_EQ(grid.edge(37), Grid::Edge (grid.node(13), grid.node(22)));
    ASSERT_EQ(grid.edge(38), Grid::Edge (grid.node(14), grid.node(23)));
    ASSERT_EQ(grid.edge(39), Grid::Edge (grid.node(15), grid.node(24)));
    ASSERT_EQ(grid.edge(40), Grid::Edge (grid.node(16), grid.node(25)));
    ASSERT_EQ(grid.edge(41), Grid::Edge (grid.node(17), grid.node(26)));
    // third slice (2D grid)
    ASSERT_EQ(grid.edge(42), Grid::Edge (grid.node(18), grid.node(19)));
    ASSERT_EQ(grid.edge(43), Grid::Edge (grid.node(19), grid.node(20)));
    ASSERT_EQ(grid.edge(44), Grid::Edge (grid.node(18), grid.node(21)));
    ASSERT_EQ(grid.edge(45), Grid::Edge (grid.node(19), grid.node(22)));
    ASSERT_EQ(grid.edge(46), Grid::Edge (grid.node(20), grid.node(23)));
    ASSERT_EQ(grid.edge(47), Grid::Edge (grid.node(21), grid.node(22)));
    ASSERT_EQ(grid.edge(48), Grid::Edge (grid.node(22), grid.node(23)));
    ASSERT_EQ(grid.edge(49), Grid::Edge (grid.node(21), grid.node(24)));
    ASSERT_EQ(grid.edge(50), Grid::Edge (grid.node(22), grid.node(25)));
    ASSERT_EQ(grid.edge(51), Grid::Edge (grid.node(23), grid.node(26)));
    ASSERT_EQ(grid.edge(52), Grid::Edge (grid.node(24), grid.node(25)));
    ASSERT_EQ(grid.edge(53), Grid::Edge (grid.node(25), grid.node(26)));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
