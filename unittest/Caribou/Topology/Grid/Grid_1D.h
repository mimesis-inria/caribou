#ifndef CARIBOU_TOPOLOGY_TEST_GRID_1D_H
#define CARIBOU_TOPOLOGY_TEST_GRID_1D_H

#include "../topology_test.h"
#include <Caribou/Topology/Grid/Grid.h>

TEST(Topology_Grid_1D, Grid1D) {

    using namespace caribou::topology;
    using Grid = Grid<1>;
    using Edge = Grid::Element;
    using WorldCoordinates = Grid::WorldCoordinates;
    using Subdivisions = Grid::Subdivisions;
    using Dimensions = Grid::Dimensions;

    Grid grid(WorldCoordinates(0.25), Subdivisions(2), Dimensions(100));

// General properties
    EXPECT_EQ(grid.number_of_nodes(), (unsigned) 3);

// Cell size test
    EXPECT_DOUBLE_EQ(100 / 2., grid.H()[0]);

// Cell numbering test
    EXPECT_EQ((Grid::CellIndex) 5, grid.cell_index_at((Grid::GridCoordinates) 5));
    EXPECT_MATRIX_EQUAL((Grid::GridCoordinates) 5, grid.cell_coordinates_at((Grid::CellIndex) 5));

// Cell positioning
    EXPECT_FALSE(grid.contains(WorldCoordinates{0.24}));
    EXPECT_TRUE(grid.contains(WorldCoordinates{0.25}));
    EXPECT_TRUE(grid.contains(WorldCoordinates{50}));
    EXPECT_TRUE(grid.contains(WorldCoordinates{100.25}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{100.26}));

    EXPECT_MATRIX_EQUAL((Grid::GridCoordinates) 0, grid.cell_coordinates_at(WorldCoordinates(0.25)));
    EXPECT_MATRIX_EQUAL((Grid::GridCoordinates) 0, grid.cell_coordinates_at(WorldCoordinates(50.24)));
    EXPECT_MATRIX_EQUAL((Grid::GridCoordinates) 1, grid.cell_coordinates_at(WorldCoordinates(50.25)));
    EXPECT_MATRIX_EQUAL((Grid::GridCoordinates) 1, grid.cell_coordinates_at(WorldCoordinates(100)));

// Cells queries
    EXPECT_LIST_EQUAL((std::list<int>{}),     grid.cells_around(WorldCoordinates(-50)));
    EXPECT_LIST_EQUAL((std::list<int>{0}),    grid.cells_around(WorldCoordinates(0.25)));
    EXPECT_LIST_EQUAL((std::list<int>{0}),    grid.cells_around(WorldCoordinates(25)));
    EXPECT_LIST_EQUAL((std::list<int>{0}),    grid.cells_around(WorldCoordinates(50.249)));
    EXPECT_LIST_EQUAL((std::list<int>{0, 1}), grid.cells_around(WorldCoordinates(50.25)));
    EXPECT_LIST_EQUAL((std::list<int>{1}),    grid.cells_around(WorldCoordinates(50.251)));
    EXPECT_LIST_EQUAL((std::list<int>{1}),    grid.cells_around(WorldCoordinates(75)));
    EXPECT_LIST_EQUAL((std::list<int>{1}),    grid.cells_around(WorldCoordinates(100.25)));

    EXPECT_LIST_EQUAL((std::list<int>{}),     grid.cells_enclosing(WorldCoordinates( -50), WorldCoordinates( 0), WorldCoordinates( 0.24)));
    EXPECT_LIST_EQUAL((std::list<int>{0}),    grid.cells_enclosing(WorldCoordinates( -50), WorldCoordinates(25), WorldCoordinates(50.24)));
    EXPECT_LIST_EQUAL((std::list<int>{0}),    grid.cells_enclosing(WorldCoordinates(0.25), WorldCoordinates(25), WorldCoordinates(50.24)));
    EXPECT_LIST_EQUAL((std::list<int>{0, 1}), grid.cells_enclosing(WorldCoordinates(50.25), WorldCoordinates(75), WorldCoordinates(100.24)));
    EXPECT_LIST_EQUAL((std::list<int>{0, 1}), grid.cells_enclosing(WorldCoordinates(25), WorldCoordinates(75)));
    // EXPECT_TRUE(list_are_equals({0,1}, grid.cells_enclosing(WorldCoordinates(0.25),  WorldCoordinates(101))));

// Node position queries
    EXPECT_MATRIX_NEAR(grid.node(0), WorldCoordinates(0.25), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(1), WorldCoordinates(50.25), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(2), WorldCoordinates(100.25), 1e-15);

// Cell Node indices
    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_cells(); ++i) {
        const auto node_indices = grid.node_indices_of((Grid::CellIndex) i);
        const Edge e(grid.node(node_indices[0]), grid.node(node_indices[1]));
        EXPECT_MATRIX_NEAR(e.nodes(), grid.cell_at((Grid::CellIndex) i).nodes(), 1e-15);
    }

}

#endif //CARIBOU_TOPOLOGY_TEST_GRID_1D_H
