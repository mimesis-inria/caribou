#include <gtest/gtest.h>
#include <Caribou/Topology/Grid/Grid.h>
#include <Caribou/Topology/HashGrid.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Quad.h>
#include <Eigen/Core>

template<template<typename, typename...> typename H, typename T, typename... Ts>
bool list_are_equals (const std::list<T> & l1, const H<T, Ts...> & l2)
{
    if (l1.size() != l2.size())
        return false;

    typename std::list<T>::const_iterator it1;
    typename H<T, Ts...>::const_iterator it2;

    for (it1 = l1.begin(), it2 = l2.begin(); it1 != l1.end(); ++it1, ++it2)
        if (*it1 != *it2)
            return false;
    return true;
}

template <typename Node>
void ASSERT_NODE_EQ(const Node & n1, const Node & n2) {
    ASSERT_FLOAT_EQ((n2-n1).norm(), 0.);
}

template <typename Element>
void ASSERT_ELEMENTS_EQ(const Element & e1, const Element & e2) {
    for (std::size_t i = 0; i < Element::NumberOfNodes; ++i) {
        ASSERT_NODE_EQ(e1.node(i), e2.node(i));
    }
}

template <typename Element1, typename Element2>
void ASSERT_ELEMENTS_EQ(const Element1 & e1, const Element2 & e2) {
    static_assert(caribou::traits<Element1>::NumberOfNodes == caribou::traits<Element2>::NumberOfNodes);
    for (std::size_t i = 0; i < Element1::NumberOfNodes; ++i) {
        ASSERT_NODE_EQ(e1.node(i), e2.node(i));
    }
}

template <typename Element, typename... Elements>
void ASSERT_ELEMENTS_EQ(const Element & e1, const Element & e2, Elements... elements) {
    ASSERT_ELEMENTS_EQ(e1, e2);
    ASSERT_ELEMENTS_EQ(e2, std::forward<Elements>(elements)...);
}

TEST(Topology, Grid1D) {
    using namespace caribou::topology;
    using Grid = Grid<1>;
    using Edge = Grid::Element;
    using WorldCoordinates = Grid::WorldCoordinates;
    using Subdivisions = Grid::Subdivisions;
    using Dimensions = Grid::Dimensions;

    Grid grid(WorldCoordinates(0.25), Subdivisions(2), Dimensions(100));

    // General properties
    ASSERT_EQ(grid.number_of_nodes(), (unsigned) 3);

    // Cell size test
    ASSERT_FLOAT_EQ(100/2., grid.H()[0]);

    // Cell numbering test
    ASSERT_EQ((Grid::CellIndex) 5, grid.cell_index_at((Grid::GridCoordinates) 5));
    ASSERT_EQ((Grid::GridCoordinates) 5, grid.grid_coordinates_at((Grid::CellIndex) 5));

    // Positioning
    ASSERT_FALSE(grid.contains(WorldCoordinates {0.24}));
    ASSERT_TRUE(grid.contains(WorldCoordinates  {0.25}));
    ASSERT_TRUE(grid.contains(WorldCoordinates  {50}));
    ASSERT_TRUE(grid.contains(WorldCoordinates  {100.25}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {100.26}));

    ASSERT_EQ((Grid::GridCoordinates) 0, grid.grid_coordinates_at(WorldCoordinates (0.25)));
    ASSERT_EQ((Grid::GridCoordinates) 0, grid.grid_coordinates_at(WorldCoordinates (50.24)));
    ASSERT_EQ((Grid::GridCoordinates) 1, grid.grid_coordinates_at(WorldCoordinates (50.25)));
    ASSERT_EQ((Grid::GridCoordinates) 1, grid.grid_coordinates_at(WorldCoordinates (100)));

    // Cells queries
    ASSERT_TRUE(list_are_equals({},    grid.cells_enclosing(WorldCoordinates(-50),   WorldCoordinates(0),  WorldCoordinates(0.24))));
    ASSERT_TRUE(list_are_equals({0},   grid.cells_enclosing(WorldCoordinates(-50),   WorldCoordinates(25), WorldCoordinates(50.24))));
    ASSERT_TRUE(list_are_equals({0},   grid.cells_enclosing(WorldCoordinates(0.25),  WorldCoordinates(25), WorldCoordinates(50.24))));
    ASSERT_TRUE(list_are_equals({0,1},   grid.cells_enclosing(WorldCoordinates(50.25), WorldCoordinates(75), WorldCoordinates(100.24))));
    ASSERT_TRUE(list_are_equals({0,1}, grid.cells_enclosing(WorldCoordinates(25),    WorldCoordinates(75))));
//    ASSERT_TRUE(list_are_equals({0,1}, grid.cells_enclosing(WorldCoordinates(0.25),  WorldCoordinates(101))));

    ASSERT_TRUE(list_are_equals({}, grid.cells_around(WorldCoordinates(-50))));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(0.25))));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(25))));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(50.249))));
    ASSERT_TRUE(list_are_equals({0,1}, grid.cells_around(WorldCoordinates(50.25))));
    ASSERT_TRUE(list_are_equals({1}, grid.cells_around(WorldCoordinates(50.251))));
    ASSERT_TRUE(list_are_equals({1}, grid.cells_around(WorldCoordinates(75))));
    ASSERT_TRUE(list_are_equals({1}, grid.cells_around(WorldCoordinates(100.25))));

    // Node queries
    ASSERT_NODE_EQ(grid.node(0), WorldCoordinates(0.25));
    ASSERT_NODE_EQ(grid.node(1), WorldCoordinates(50.25));
    ASSERT_NODE_EQ(grid.node(2), WorldCoordinates(100.25));

    // Edge queries
    ASSERT_EQ(grid.number_of_edges(), (unsigned) 2);
    ASSERT_ELEMENTS_EQ(grid.edge(0), Grid::Edge(WorldCoordinates {0.25}, WorldCoordinates{50.25}));
    ASSERT_ELEMENTS_EQ(grid.edge(1), Grid::Edge(WorldCoordinates{50.25}, WorldCoordinates{100.25}));

    // Node indices
    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_cells(); ++i) {
        const auto node_indices = grid.node_indices_of((Grid::CellIndex) i);
        const Edge e (grid.node(node_indices[0]), grid.node(node_indices[1]));
        ASSERT_ELEMENTS_EQ(e, grid.cell_at((Grid::CellIndex) i));
    }
}

TEST(Topology, Grid2D) {
    using namespace caribou::topology;
    using Grid = Grid<2>;

    using WorldCoordinates = Grid::WorldCoordinates;
    using Subdivisions = Grid::Subdivisions;
    using Dimensions = Grid::Dimensions;
    using GridCoordinates = Grid::GridCoordinates;
    using CellIndex = Grid::CellIndex;
    using Edge = Grid::Edge;
    using Quad = caribou::geometry::Quad<2>;

    Grid grid(WorldCoordinates {0.25, 0.5}, Subdivisions {2, 2}, Dimensions {100, 100});

    // General properties
    ASSERT_EQ(grid.number_of_nodes(), (unsigned) 9);

    // Cell numbering test
    ASSERT_EQ((CellIndex) 3, grid.cell_index_at(GridCoordinates (1,1)));
    ASSERT_EQ(GridCoordinates (1,1), grid.grid_coordinates_at((CellIndex) 3));

    // Positioning
    ASSERT_FALSE(grid.contains(WorldCoordinates {0.24,0.5}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {0.25,0.49}));
    ASSERT_TRUE (grid.contains(WorldCoordinates {0.25,0.5}));
    ASSERT_TRUE (grid.contains(WorldCoordinates {50,50}));
    ASSERT_TRUE (grid.contains(WorldCoordinates {100.25,100.5}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {100.25,100.51}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {100.26,100.5}));

    // Cells queries
    ASSERT_TRUE(list_are_equals({0,2}, grid.cells_enclosing(
            WorldCoordinates ({-50, 50}),
            WorldCoordinates ({25, 25}),
            WorldCoordinates ({25, 75})
    )));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_enclosing(
            WorldCoordinates ({0.25, 0.5}),
            WorldCoordinates ({25, 25}),
            WorldCoordinates ({50.24, 50.49})
            )));
    ASSERT_TRUE(list_are_equals({0, 1}, grid.cells_enclosing(
            WorldCoordinates ({50.25, 0.5}),
            WorldCoordinates ({75, 25}),
            WorldCoordinates ({100.24, 50.49})
    )));
    ASSERT_TRUE(list_are_equals({0,1,2,3}, grid.cells_enclosing(
            WorldCoordinates ({0.25, 0.5}),
            WorldCoordinates ({75.25, 75.50})
    )));

    ASSERT_TRUE(list_are_equals({}, grid.cells_around(WorldCoordinates(-50, -50))));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(0.25, 0.50))));
    ASSERT_TRUE(list_are_equals({0, 1}, grid.cells_around(WorldCoordinates(50.25, 0.50))));
    ASSERT_TRUE(list_are_equals({1}, grid.cells_around(WorldCoordinates(100.25, 0.50))));
    ASSERT_TRUE(list_are_equals({0, 2}, grid.cells_around(WorldCoordinates(0.25, 50.50))));
    ASSERT_TRUE(list_are_equals({0, 2, 1, 3}, grid.cells_around(WorldCoordinates(50.25, 50.50))));
    ASSERT_TRUE(list_are_equals({1, 3}, grid.cells_around(WorldCoordinates(100.25, 50.50))));
    ASSERT_TRUE(list_are_equals({2}, grid.cells_around(WorldCoordinates(0.25, 100.50))));
    ASSERT_TRUE(list_are_equals({2, 3}, grid.cells_around(WorldCoordinates(50.25, 100.50))));
    ASSERT_TRUE(list_are_equals({3}, grid.cells_around(WorldCoordinates(100.25, 100.50))));

    // Node queries
    ASSERT_NODE_EQ(grid.node(0), WorldCoordinates ({0.25, 0.5}));
    ASSERT_NODE_EQ(grid.node(1), WorldCoordinates ({50.25, 0.5}));
    ASSERT_NODE_EQ(grid.node(2), WorldCoordinates ({100.25, 0.5}));
    ASSERT_NODE_EQ(grid.node(3), WorldCoordinates ({0.25, 50.5}));
    ASSERT_NODE_EQ(grid.node(4), WorldCoordinates ({50.25, 50.5}));
    ASSERT_NODE_EQ(grid.node(5), WorldCoordinates ({100.25, 50.5}));
    ASSERT_NODE_EQ(grid.node(6), WorldCoordinates ({0.25, 100.5}));
    ASSERT_NODE_EQ(grid.node(7), WorldCoordinates ({50.25, 100.5}));
    ASSERT_NODE_EQ(grid.node(8), WorldCoordinates ({100.25, 100.5}));

    // Edge queries
    ASSERT_EQ(grid.number_of_edges(), (unsigned) 12);
    ASSERT_ELEMENTS_EQ(grid.edge(0),  Edge (grid.node(0), grid.node(1)));
    ASSERT_ELEMENTS_EQ(grid.edge(1),  Edge (grid.node(1), grid.node(2)));
    ASSERT_ELEMENTS_EQ(grid.edge(2),  Edge (grid.node(0), grid.node(3)));
    ASSERT_ELEMENTS_EQ(grid.edge(3),  Edge (grid.node(1), grid.node(4)));
    ASSERT_ELEMENTS_EQ(grid.edge(4),  Edge (grid.node(2), grid.node(5)));
    ASSERT_ELEMENTS_EQ(grid.edge(5),  Edge (grid.node(3), grid.node(4)));
    ASSERT_ELEMENTS_EQ(grid.edge(6),  Edge (grid.node(4), grid.node(5)));
    ASSERT_ELEMENTS_EQ(grid.edge(7),  Edge (grid.node(3), grid.node(6)));
    ASSERT_ELEMENTS_EQ(grid.edge(8),  Edge (grid.node(4), grid.node(7)));
    ASSERT_ELEMENTS_EQ(grid.edge(9),  Edge (grid.node(5), grid.node(8)));
    ASSERT_ELEMENTS_EQ(grid.edge(10), Edge (grid.node(6), grid.node(7)));
    ASSERT_ELEMENTS_EQ(grid.edge(11), Edge (grid.node(7), grid.node(8)));

    // Node indices
    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_cells(); ++i) {
        const auto node_indices = grid.node_indices_of((Grid::CellIndex) i);
        const Quad q (grid.node(node_indices[0]), grid.node(node_indices[1]), grid.node(node_indices[2]), grid.node(node_indices[3]));
        ASSERT_ELEMENTS_EQ(q, grid.cell_at((Grid::CellIndex) i));
    }
}

TEST(Topology, Grid3D) {
    using namespace caribou::topology;
    using Grid = Grid<3>;

    using WorldCoordinates = Grid::WorldCoordinates;
    using Subdivisions = Grid::Subdivisions;
    using Dimensions = Grid::Dimensions;
    using GridCoordinates = Grid::GridCoordinates;
    using CellIndex = Grid::CellIndex;
    using Edge = Grid::Edge;
    using Hexa = caribou::geometry::Hexahedron<caribou::geometry::interpolation::Hexahedron8>;

    Grid grid(WorldCoordinates {0.25, 0.5, 0.75}, Subdivisions {2, 2, 2}, Dimensions {100, 100, 100});

    // General properties
    ASSERT_EQ(grid.number_of_nodes(), (unsigned) 27);

    // Cell numbering test
    ASSERT_EQ((CellIndex) 5, grid.cell_index_at(GridCoordinates (1,0,1)));
    ASSERT_EQ(GridCoordinates (1,0,1), grid.grid_coordinates_at((CellIndex) 5));

    // Positioning
    ASSERT_FALSE(grid.contains(WorldCoordinates {0.24,0.50, 0.75}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {0.25,0.49,0.75}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {0.25,0.50,0.74}));
    ASSERT_TRUE (grid.contains(WorldCoordinates {0.25,0.50,0.75}));
    ASSERT_TRUE (grid.contains(WorldCoordinates {50,50,50}));
    ASSERT_TRUE (grid.contains(WorldCoordinates {100.25,100.50,100.75}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {100.26,100.50,100.75}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {100.25,100.51,100.75}));
    ASSERT_FALSE(grid.contains(WorldCoordinates {100.25,100.50,100.76}));

    // Cells queries
    ASSERT_TRUE(list_are_equals({0,2,4,6}, grid.cells_enclosing(
            WorldCoordinates ({-50, 50, 50}),
            WorldCoordinates ({25, 25, 25}),
            WorldCoordinates ({25, 75, 75})
    )));

    ASSERT_TRUE(list_are_equals({0}, grid.cells_enclosing(
            WorldCoordinates ({0.25, 0.5, 0.75}),
            WorldCoordinates ({25, 25, 25}),
            WorldCoordinates ({50.24, 50.49, 50.74})
    )));
    ASSERT_TRUE(list_are_equals({0,1,2,3,4,5,6,7}, grid.cells_enclosing(
            WorldCoordinates ({50.25, 50.5, 50.75}),
            WorldCoordinates ({75, 75, 75}),
            WorldCoordinates ({100.24, 100.49, 100.74})
    )));
    ASSERT_TRUE(list_are_equals({0,1,2,3,4,5,6,7}, grid.cells_enclosing(
            WorldCoordinates ({0.25, 0.5, 0.75}),
            WorldCoordinates ({100.24, 100.49, 100.74})
    )));

    // Cells around nodes
    ASSERT_TRUE(list_are_equals({}, grid.cells_around(WorldCoordinates(-50, -50, -50))));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(0.25, 0.50, 0.75))));
    ASSERT_TRUE(list_are_equals({0, 1}, grid.cells_around(WorldCoordinates(50.25, 0.50, 0.75))));
    ASSERT_TRUE(list_are_equals({1}, grid.cells_around(WorldCoordinates(100.25, 0.50, 0.75))));
    ASSERT_TRUE(list_are_equals({0, 4, 2, 6, 1, 5, 3, 7}, grid.cells_around(WorldCoordinates(50.25, 50.50, 50.75))));

    // Cells around faces
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(0.25, 25.50, 25.75))));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(25.25, 25.50, 0.75))));
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(25.25, 0.50, 25.75))));
    ASSERT_TRUE(list_are_equals({0, 1}, grid.cells_around(WorldCoordinates(50.25, 25.50, 25.75))));
    ASSERT_TRUE(list_are_equals({0, 2}, grid.cells_around(WorldCoordinates(25.25, 50.50, 25.75))));
    ASSERT_TRUE(list_are_equals({0, 4}, grid.cells_around(WorldCoordinates(25.25, 25.50, 50.75))));
    ASSERT_TRUE(list_are_equals({7}, grid.cells_around(WorldCoordinates(75.25, 100.50, 75.75))));
    ASSERT_TRUE(list_are_equals({7}, grid.cells_around(WorldCoordinates(75.25, 75.50, 100.75))));
    ASSERT_TRUE(list_are_equals({7}, grid.cells_around(WorldCoordinates(100.25, 75.50, 75.75))));

    // Cells around edges
    ASSERT_TRUE(list_are_equals({0}, grid.cells_around(WorldCoordinates(25.25, 0.50, 0.75))));
    ASSERT_TRUE(list_are_equals({7}, grid.cells_around(WorldCoordinates(75.25, 100.50, 100.75))));
    ASSERT_TRUE(list_are_equals({0, 4, 2, 6}, grid.cells_around(WorldCoordinates(25.25, 50.50, 50.75))));
    ASSERT_TRUE(list_are_equals({2, 6, 3, 7}, grid.cells_around(WorldCoordinates(50.25, 75.50, 50.75))));
    ASSERT_TRUE(list_are_equals({0, 2, 1, 3}, grid.cells_around(WorldCoordinates(50.25, 50.50, 25.75))));

    // Node queries
    ASSERT_NODE_EQ(grid.node(0), WorldCoordinates ({0.25, 0.5, 0.75}));
    ASSERT_NODE_EQ(grid.node(2), WorldCoordinates ({100.25, 0.5, 0.75}));
    ASSERT_NODE_EQ(grid.node(6), WorldCoordinates ({0.25, 100.5, 0.75}));
    ASSERT_NODE_EQ(grid.node(8), WorldCoordinates ({100.25, 100.5, 0.75}));

    ASSERT_NODE_EQ(grid.node(18), WorldCoordinates ({0.25, 0.5, 100.75}));
    ASSERT_NODE_EQ(grid.node(20), WorldCoordinates ({100.25, 0.5, 100.75}));
    ASSERT_NODE_EQ(grid.node(24), WorldCoordinates ({0.25, 100.5, 100.75}));
    ASSERT_NODE_EQ(grid.node(26), WorldCoordinates ({100.25, 100.5, 100.75}));

    // Edge queries
    ASSERT_EQ(grid.number_of_edges(), (unsigned) 3*12+2*9);
    // First slice (2D grid)
    ASSERT_ELEMENTS_EQ(grid.edge(0),  Edge (grid.node(0), grid.node(1)));
    ASSERT_ELEMENTS_EQ(grid.edge(1),  Edge (grid.node(1), grid.node(2)));
    ASSERT_ELEMENTS_EQ(grid.edge(2),  Edge (grid.node(0), grid.node(3)));
    ASSERT_ELEMENTS_EQ(grid.edge(3),  Edge (grid.node(1), grid.node(4)));
    ASSERT_ELEMENTS_EQ(grid.edge(4),  Edge (grid.node(2), grid.node(5)));
    ASSERT_ELEMENTS_EQ(grid.edge(5),  Edge (grid.node(3), grid.node(4)));
    ASSERT_ELEMENTS_EQ(grid.edge(6),  Edge (grid.node(4), grid.node(5)));
    ASSERT_ELEMENTS_EQ(grid.edge(7),  Edge (grid.node(3), grid.node(6)));
    ASSERT_ELEMENTS_EQ(grid.edge(8),  Edge (grid.node(4), grid.node(7)));
    ASSERT_ELEMENTS_EQ(grid.edge(9),  Edge (grid.node(5), grid.node(8)));
    ASSERT_ELEMENTS_EQ(grid.edge(10), Edge (grid.node(6), grid.node(7)));
    ASSERT_ELEMENTS_EQ(grid.edge(11), Edge (grid.node(7), grid.node(8)));
    // Between the slice 1 and slice 2
    ASSERT_ELEMENTS_EQ(grid.edge(12), Edge (grid.node(0), grid.node(9)));
    ASSERT_ELEMENTS_EQ(grid.edge(13), Edge (grid.node(1), grid.node(10)));
    ASSERT_ELEMENTS_EQ(grid.edge(14), Edge (grid.node(2), grid.node(11)));
    ASSERT_ELEMENTS_EQ(grid.edge(15), Edge (grid.node(3), grid.node(12)));
    ASSERT_ELEMENTS_EQ(grid.edge(16), Edge (grid.node(4), grid.node(13)));
    ASSERT_ELEMENTS_EQ(grid.edge(17), Edge (grid.node(5), grid.node(14)));
    ASSERT_ELEMENTS_EQ(grid.edge(18), Edge (grid.node(6), grid.node(15)));
    ASSERT_ELEMENTS_EQ(grid.edge(19), Edge (grid.node(7), grid.node(16)));
    ASSERT_ELEMENTS_EQ(grid.edge(20), Edge (grid.node(8), grid.node(17)));
    // second slice (2D grid)
    ASSERT_ELEMENTS_EQ(grid.edge(21), Edge (grid.node(9), grid.node(10)));
    ASSERT_ELEMENTS_EQ(grid.edge(22), Edge (grid.node(10), grid.node(11)));
    ASSERT_ELEMENTS_EQ(grid.edge(23), Edge (grid.node(9), grid.node(12)));
    ASSERT_ELEMENTS_EQ(grid.edge(24), Edge (grid.node(10), grid.node(13)));
    ASSERT_ELEMENTS_EQ(grid.edge(25), Edge (grid.node(11), grid.node(14)));
    ASSERT_ELEMENTS_EQ(grid.edge(26), Edge (grid.node(12), grid.node(13)));
    ASSERT_ELEMENTS_EQ(grid.edge(27), Edge (grid.node(13), grid.node(14)));
    ASSERT_ELEMENTS_EQ(grid.edge(28), Edge (grid.node(12), grid.node(15)));
    ASSERT_ELEMENTS_EQ(grid.edge(29), Edge (grid.node(13), grid.node(16)));
    ASSERT_ELEMENTS_EQ(grid.edge(30), Edge (grid.node(14), grid.node(17)));
    ASSERT_ELEMENTS_EQ(grid.edge(31), Edge (grid.node(15), grid.node(16)));
    ASSERT_ELEMENTS_EQ(grid.edge(32), Edge (grid.node(16), grid.node(17)));
    // Between the slice 2 and slice 3
    ASSERT_ELEMENTS_EQ(grid.edge(33), Edge (grid.node(9), grid.node(18)));
    ASSERT_ELEMENTS_EQ(grid.edge(34), Edge (grid.node(10), grid.node(19)));
    ASSERT_ELEMENTS_EQ(grid.edge(35), Edge (grid.node(11), grid.node(20)));
    ASSERT_ELEMENTS_EQ(grid.edge(36), Edge (grid.node(12), grid.node(21)));
    ASSERT_ELEMENTS_EQ(grid.edge(37), Edge (grid.node(13), grid.node(22)));
    ASSERT_ELEMENTS_EQ(grid.edge(38), Edge (grid.node(14), grid.node(23)));
    ASSERT_ELEMENTS_EQ(grid.edge(39), Edge (grid.node(15), grid.node(24)));
    ASSERT_ELEMENTS_EQ(grid.edge(40), Edge (grid.node(16), grid.node(25)));
    ASSERT_ELEMENTS_EQ(grid.edge(41), Edge (grid.node(17), grid.node(26)));
    // third slice (2D grid)
    ASSERT_ELEMENTS_EQ(grid.edge(42), Edge (grid.node(18), grid.node(19)));
    ASSERT_ELEMENTS_EQ(grid.edge(43), Edge (grid.node(19), grid.node(20)));
    ASSERT_ELEMENTS_EQ(grid.edge(44), Edge (grid.node(18), grid.node(21)));
    ASSERT_ELEMENTS_EQ(grid.edge(45), Edge (grid.node(19), grid.node(22)));
    ASSERT_ELEMENTS_EQ(grid.edge(46), Edge (grid.node(20), grid.node(23)));
    ASSERT_ELEMENTS_EQ(grid.edge(47), Edge (grid.node(21), grid.node(22)));
    ASSERT_ELEMENTS_EQ(grid.edge(48), Edge (grid.node(22), grid.node(23)));
    ASSERT_ELEMENTS_EQ(grid.edge(49), Edge (grid.node(21), grid.node(24)));
    ASSERT_ELEMENTS_EQ(grid.edge(50), Edge (grid.node(22), grid.node(25)));
    ASSERT_ELEMENTS_EQ(grid.edge(51), Edge (grid.node(23), grid.node(26)));
    ASSERT_ELEMENTS_EQ(grid.edge(52), Edge (grid.node(24), grid.node(25)));
    ASSERT_ELEMENTS_EQ(grid.edge(53), Edge (grid.node(25), grid.node(26)));

    // Node indices
    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_cells(); ++i) {
        const auto node_indices = grid.node_indices_of((Grid::CellIndex) i);
        const Hexa h (
            grid.node(node_indices[0]), grid.node(node_indices[1]), grid.node(node_indices[2]), grid.node(node_indices[3]),
            grid.node(node_indices[4]), grid.node(node_indices[5]), grid.node(node_indices[6]), grid.node(node_indices[7]));
        ASSERT_ELEMENTS_EQ(h, grid.cell_at((Grid::CellIndex) i));
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
