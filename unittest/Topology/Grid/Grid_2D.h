#ifndef CARIBOU_TOPOLOGY_TEST_GRID_2D_H
#define CARIBOU_TOPOLOGY_TEST_GRID_2D_H

#include "../topology_test.h"
#include <Caribou/Topology/Grid/Grid.h>
#include <Caribou/Geometry/Quad.h>

TEST(Topology_Grid_2D, Grid2D) {
    using namespace caribou::topology;
    using Grid = Grid<2>;

    using WorldCoordinates = Grid::WorldCoordinates;
    using Subdivisions = Grid::Subdivisions;
    using Dimensions = Grid::Dimensions;
    using GridCoordinates = Grid::GridCoordinates;
    using CellIndex = Grid::CellIndex;
    using Quad = caribou::geometry::Quad<2>;

    Grid grid(WorldCoordinates{0.25, 0.5}, Subdivisions{2, 2}, Dimensions{100, 100});

// General properties
    EXPECT_EQ(grid.number_of_nodes(), (unsigned) 9);

// Cell numbering test
    EXPECT_EQ((CellIndex) 3, grid.cell_index_at({1, 1}));
    EXPECT_MATRIX_EQUAL(GridCoordinates(1, 1), grid.cell_coordinates_at((CellIndex) 3));

// Cell positioning
    EXPECT_FALSE(grid.contains(WorldCoordinates{0.24, 0.5}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{0.25, 0.49}));
    EXPECT_TRUE (grid.contains(WorldCoordinates{0.25, 0.5}));
    EXPECT_TRUE (grid.contains(WorldCoordinates{50, 50}));
    EXPECT_TRUE (grid.contains(WorldCoordinates{100.25, 100.5}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{100.25, 100.51}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{100.26, 100.5}));

// Cells queries
    EXPECT_LIST_EQUAL((std::list({0, 2})), grid.cells_enclosing(
        WorldCoordinates({-50, 50}),
        WorldCoordinates({25, 25}),
        WorldCoordinates({25, 75})
    ));
    EXPECT_LIST_EQUAL((std::list({0})), grid.cells_enclosing(
        WorldCoordinates({0.25, 0.5}),
        WorldCoordinates({25, 25}),
        WorldCoordinates({50.24, 50.49})
    ));
    EXPECT_LIST_EQUAL((std::list{0, 1}), grid.cells_enclosing(
        WorldCoordinates({50.25, 0.5}),
        WorldCoordinates({75, 25}),
        WorldCoordinates({100.24, 50.49})
    ));
    EXPECT_LIST_EQUAL((std::list{0, 1, 2, 3}), grid.cells_enclosing(
        WorldCoordinates({0.25, 0.5}),
        WorldCoordinates({75.25, 75.50})
    ));

    EXPECT_LIST_EQUAL((std::list<int>{}),           grid.cells_around(WorldCoordinates(-50, -50)));
    EXPECT_LIST_EQUAL((std::list<int>{0}),          grid.cells_around(WorldCoordinates(0.25, 0.50)));
    EXPECT_LIST_EQUAL((std::list<int>{0, 1}),       grid.cells_around(WorldCoordinates(50.25, 0.50)));
    EXPECT_LIST_EQUAL((std::list<int>{1}),          grid.cells_around(WorldCoordinates(100.25, 0.50)));
    EXPECT_LIST_EQUAL((std::list<int>{0, 2}),       grid.cells_around(WorldCoordinates(0.25, 50.50)));
    EXPECT_LIST_EQUAL((std::list<int>{0, 2, 1, 3}), grid.cells_around(WorldCoordinates(50.25, 50.50)));
    EXPECT_LIST_EQUAL((std::list<int>{1, 3}),       grid.cells_around(WorldCoordinates(100.25, 50.50)));
    EXPECT_LIST_EQUAL((std::list<int>{2}),          grid.cells_around(WorldCoordinates(0.25, 100.50)));
    EXPECT_LIST_EQUAL((std::list<int>{2, 3}),       grid.cells_around(WorldCoordinates(50.25, 100.50)));
    EXPECT_LIST_EQUAL((std::list<int>{3}),          grid.cells_around(WorldCoordinates(100.25, 100.50)));

// Node position queries
    EXPECT_MATRIX_NEAR(grid.node(0), WorldCoordinates({0.25, 0.5}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(1), WorldCoordinates({50.25, 0.5}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(2), WorldCoordinates({100.25, 0.5}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(3), WorldCoordinates({0.25, 50.5}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(4), WorldCoordinates({50.25, 50.5}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(5), WorldCoordinates({100.25, 50.5}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(6), WorldCoordinates({0.25, 100.5}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(7), WorldCoordinates({50.25, 100.5}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(8), WorldCoordinates({100.25, 100.5}), 1e-15);

// Node indexing
    for (Grid::NodeIndex index = 0; index < (signed) grid.number_of_nodes(); ++index) {
        EXPECT_EQ(index, grid.node_index_at(grid.node_coordinates_at(index)));
    }

// Edge queries
    EXPECT_EQ(grid.number_of_edges(), (unsigned) 12);
    EXPECT_EQ(grid.edge(0), Grid::EdgeNodes({{0, 1}}));
    EXPECT_EQ(grid.edge(1), Grid::EdgeNodes({{1, 2}}));
    EXPECT_EQ(grid.edge(2), Grid::EdgeNodes({{0, 3}}));
    EXPECT_EQ(grid.edge(3), Grid::EdgeNodes({{1, 4}}));
    EXPECT_EQ(grid.edge(4), Grid::EdgeNodes({{2, 5}}));
    EXPECT_EQ(grid.edge(5), Grid::EdgeNodes({{3, 4}}));
    EXPECT_EQ(grid.edge(6), Grid::EdgeNodes({{4, 5}}));
    EXPECT_EQ(grid.edge(7), Grid::EdgeNodes({{3, 6}}));
    EXPECT_EQ(grid.edge(8), Grid::EdgeNodes({{4, 7}}));
    EXPECT_EQ(grid.edge(9), Grid::EdgeNodes({{5, 8}}));
    EXPECT_EQ(grid.edge(10), Grid::EdgeNodes({{6, 7}}));
    EXPECT_EQ(grid.edge(11), Grid::EdgeNodes({{7, 8}}));

// Cell Node indices
    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_cells(); ++i) {
        const auto node_indices = grid.node_indices_of((Grid::CellIndex) i);
        const Quad q(grid.node(node_indices[0]), grid.node(node_indices[1]), grid.node(node_indices[2]),
                     grid.node(node_indices[3]));
        EXPECT_MATRIX_NEAR(q.nodes(), grid.cell_at((Grid::CellIndex) i).nodes(), 1e-15);
    }
}

#endif //CARIBOU_TOPOLOGY_TEST_GRID_2D_H
