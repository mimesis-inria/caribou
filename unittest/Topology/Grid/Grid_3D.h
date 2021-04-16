#ifndef CARIBOU_TOPOLOGY_TEST_GRID_3D_H
#define CARIBOU_TOPOLOGY_TEST_GRID_3D_H

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Topology/Grid/Grid.h>
#include <vector>
#include <chrono>

#define BEGIN_CLOCK ;std::chrono::steady_clock::time_point __time_point_begin;
#define TICK ;__time_point_begin = std::chrono::steady_clock::now();
#define TOCK (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - __time_point_begin).count())

TEST(Topology_Grid_3D, Grid3D) {
    using namespace caribou::topology;
    using Grid = Grid<3>;

    using WorldCoordinates = Grid::WorldCoordinates;
    using Subdivisions = Grid::Subdivisions;
    using Dimensions = Grid::Dimensions;
    using GridCoordinates = Grid::GridCoordinates;
    using CellIndex = Grid::CellIndex;
    using Hexa = caribou::geometry::Hexahedron<caribou::Linear>;

    Grid grid(WorldCoordinates{0.25, 0.5, 0.75}, Subdivisions{2, 2, 2}, Dimensions{100, 100, 100});

    // General properties
    EXPECT_EQ(grid.number_of_nodes(), (unsigned) 27);

    // Cell numbering test
    EXPECT_EQ((CellIndex) 5, grid.cell_index_at({1, 0, 1}));
    EXPECT_MATRIX_EQUAL(GridCoordinates(1, 0, 1), grid.cell_coordinates_at((CellIndex) 5));

    // Cell positioning
    EXPECT_FALSE(grid.contains(WorldCoordinates{ 0.24,    0.50,   0.75}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{ 0.25,    0.49,   0.75}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{ 0.25,    0.50,   0.74}));
    EXPECT_TRUE (grid.contains(WorldCoordinates{ 0.25,    0.50,   0.75}));
    EXPECT_TRUE (grid.contains(WorldCoordinates{ 50.00,  50.00,  50.00}));
    EXPECT_TRUE (grid.contains(WorldCoordinates{100.25, 100.50, 100.75}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{100.26, 100.50, 100.75}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{100.25, 100.51, 100.75}));
    EXPECT_FALSE(grid.contains(WorldCoordinates{100.25, 100.50, 100.76}));

    // Cells queries
    EXPECT_LIST_EQUAL((std::list{0, 2, 4, 6}), grid.cells_enclosing(
        WorldCoordinates({-50, 50, 50}),
        WorldCoordinates({25, 25, 25}),
        WorldCoordinates({25, 75, 75})
    ));

    EXPECT_LIST_EQUAL((std::list{0}), grid.cells_enclosing(
        WorldCoordinates({0.25, 0.5, 0.75}),
        WorldCoordinates({25, 25, 25}),
        WorldCoordinates({50.24, 50.49, 50.74})
    ));
    EXPECT_LIST_EQUAL((std::list{0, 1, 2, 3, 4, 5, 6, 7}), grid.cells_enclosing(
        WorldCoordinates({50.25, 50.5, 50.75}),
        WorldCoordinates({75, 75, 75}),
        WorldCoordinates({100.24, 100.49, 100.74})
    ));
    EXPECT_LIST_EQUAL((std::list{0, 1, 2, 3, 4, 5, 6, 7}), grid.cells_enclosing(
        WorldCoordinates({0.25, 0.5, 0.75}),
        WorldCoordinates({100.24, 100.49, 100.74})
    ));

    // Cells around nodes
    EXPECT_LIST_EQUAL((std::list<int>{}), grid.cells_around(WorldCoordinates(-50, -50, -50)));
    EXPECT_LIST_EQUAL((std::list{0}),    grid.cells_around(WorldCoordinates(  0.25, 0.50, 0.75)));
    EXPECT_LIST_EQUAL((std::list{0, 1}), grid.cells_around(WorldCoordinates( 50.25, 0.50, 0.75)));
    EXPECT_LIST_EQUAL((std::list{1}),    grid.cells_around(WorldCoordinates(100.25, 0.50, 0.75)));
    EXPECT_LIST_EQUAL((std::list{0, 4, 2, 6, 1, 5, 3, 7}), grid.cells_around(WorldCoordinates(50.25, 50.50, 50.75)));

    // Cells around faces
    auto eps = std::numeric_limits<WorldCoordinates::Scalar>::min();
    EXPECT_LIST_EQUAL((std::list{0}), grid.cells_around(WorldCoordinates(0.25-eps, 25.50, 25.75)));
    EXPECT_LIST_EQUAL((std::list{0}), grid.cells_around(WorldCoordinates(25.25, 25.50, 0.75)));
    EXPECT_LIST_EQUAL((std::list{0}), grid.cells_around(WorldCoordinates(25.25, 0.50, 25.75)));
    EXPECT_LIST_EQUAL((std::list{0, 1}), grid.cells_around(WorldCoordinates(50.25, 25.50, 25.75)));
    EXPECT_LIST_EQUAL((std::list{0, 2}), grid.cells_around(WorldCoordinates(25.25, 50.50, 25.75)));
    EXPECT_LIST_EQUAL((std::list{0, 4}), grid.cells_around(WorldCoordinates(25.25, 25.50, 50.75)));
    EXPECT_LIST_EQUAL((std::list{7}), grid.cells_around(WorldCoordinates(75.25, 100.50, 75.75)));
    EXPECT_LIST_EQUAL((std::list{7}), grid.cells_around(WorldCoordinates(75.25, 75.50, 100.75)));
    EXPECT_LIST_EQUAL((std::list{7}), grid.cells_around(WorldCoordinates(100.25, 75.50, 75.75)));

    // Cells around edges
    EXPECT_LIST_EQUAL((std::list{0}), grid.cells_around(WorldCoordinates(25.25, 0.50, 0.75)));
    EXPECT_LIST_EQUAL((std::list{7}), grid.cells_around(WorldCoordinates(75.25, 100.50, 100.75)));
    EXPECT_LIST_EQUAL((std::list{0, 4, 2, 6}), grid.cells_around(WorldCoordinates(25.25, 50.50, 50.75)));
    EXPECT_LIST_EQUAL((std::list{2, 6, 3, 7}), grid.cells_around(WorldCoordinates(50.25, 75.50, 50.75)));
    EXPECT_LIST_EQUAL((std::list{0, 2, 1, 3}), grid.cells_around(WorldCoordinates(50.25, 50.50, 25.75)));

    // Node position queries
    EXPECT_MATRIX_NEAR(grid.node(0), WorldCoordinates({  0.25,   0.5, 0.75}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(2), WorldCoordinates({100.25,   0.5, 0.75}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(6), WorldCoordinates({  0.25, 100.5, 0.75}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(8), WorldCoordinates({100.25, 100.5, 0.75}), 1e-15);

    EXPECT_MATRIX_NEAR(grid.node(18), WorldCoordinates({  0.25,   0.5, 100.75}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(20), WorldCoordinates({100.25,   0.5, 100.75}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(24), WorldCoordinates({  0.25, 100.5, 100.75}), 1e-15);
    EXPECT_MATRIX_NEAR(grid.node(26), WorldCoordinates({100.25, 100.5, 100.75}), 1e-15);

    // Node indexing
    for (Grid::NodeIndex index = 0; index < (signed) grid.number_of_nodes(); ++index) {
        EXPECT_EQ(index, grid.node_index_at(grid.node_coordinates_at(index)));
    }

    // Edge queries
    EXPECT_EQ(grid.number_of_edges(), (unsigned) 3 * 12 + 2 * 9);
    // First slice (2D grid)
    EXPECT_EQ(grid.edge(0), Grid::EdgeNodes({{0, 1}}));
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
    // Between the slice 1 and slice 2
    EXPECT_EQ(grid.edge(12), Grid::EdgeNodes({{0, 9}}));
    EXPECT_EQ(grid.edge(13), Grid::EdgeNodes({{1, 10}}));
    EXPECT_EQ(grid.edge(14), Grid::EdgeNodes({{2, 11}}));
    EXPECT_EQ(grid.edge(15), Grid::EdgeNodes({{3, 12}}));
    EXPECT_EQ(grid.edge(16), Grid::EdgeNodes({{4, 13}}));
    EXPECT_EQ(grid.edge(17), Grid::EdgeNodes({{5, 14}}));
    EXPECT_EQ(grid.edge(18), Grid::EdgeNodes({{6, 15}}));
    EXPECT_EQ(grid.edge(19), Grid::EdgeNodes({{7, 16}}));
    EXPECT_EQ(grid.edge(20), Grid::EdgeNodes({{8, 17}}));
    // second slice (2D grid)
    EXPECT_EQ(grid.edge(21), Grid::EdgeNodes({{9, 10}}));
    EXPECT_EQ(grid.edge(22), Grid::EdgeNodes({{10, 11}}));
    EXPECT_EQ(grid.edge(23), Grid::EdgeNodes({{9, 12}}));
    EXPECT_EQ(grid.edge(24), Grid::EdgeNodes({{10, 13}}));
    EXPECT_EQ(grid.edge(25), Grid::EdgeNodes({{11, 14}}));
    EXPECT_EQ(grid.edge(26), Grid::EdgeNodes({{12, 13}}));
    EXPECT_EQ(grid.edge(27), Grid::EdgeNodes({{13, 14}}));
    EXPECT_EQ(grid.edge(28), Grid::EdgeNodes({{12, 15}}));
    EXPECT_EQ(grid.edge(29), Grid::EdgeNodes({{13, 16}}));
    EXPECT_EQ(grid.edge(30), Grid::EdgeNodes({{14, 17}}));
    EXPECT_EQ(grid.edge(31), Grid::EdgeNodes({{15, 16}}));
    EXPECT_EQ(grid.edge(32), Grid::EdgeNodes({{16, 17}}));
    // Between the slice 2 and slice 3
    EXPECT_EQ(grid.edge(33), Grid::EdgeNodes({{9, 18}}));
    EXPECT_EQ(grid.edge(34), Grid::EdgeNodes({{10, 19}}));
    EXPECT_EQ(grid.edge(35), Grid::EdgeNodes({{11, 20}}));
    EXPECT_EQ(grid.edge(36), Grid::EdgeNodes({{12, 21}}));
    EXPECT_EQ(grid.edge(37), Grid::EdgeNodes({{13, 22}}));
    EXPECT_EQ(grid.edge(38), Grid::EdgeNodes({{14, 23}}));
    EXPECT_EQ(grid.edge(39), Grid::EdgeNodes({{15, 24}}));
    EXPECT_EQ(grid.edge(40), Grid::EdgeNodes({{16, 25}}));
    EXPECT_EQ(grid.edge(41), Grid::EdgeNodes({{17, 26}}));
    // third slice (2D grid)
    EXPECT_EQ(grid.edge(42), Grid::EdgeNodes({{18, 19}}));
    EXPECT_EQ(grid.edge(43), Grid::EdgeNodes({{19, 20}}));
    EXPECT_EQ(grid.edge(44), Grid::EdgeNodes({{18, 21}}));
    EXPECT_EQ(grid.edge(45), Grid::EdgeNodes({{19, 22}}));
    EXPECT_EQ(grid.edge(46), Grid::EdgeNodes({{20, 23}}));
    EXPECT_EQ(grid.edge(47), Grid::EdgeNodes({{21, 22}}));
    EXPECT_EQ(grid.edge(48), Grid::EdgeNodes({{22, 23}}));
    EXPECT_EQ(grid.edge(49), Grid::EdgeNodes({{21, 24}}));
    EXPECT_EQ(grid.edge(50), Grid::EdgeNodes({{22, 25}}));
    EXPECT_EQ(grid.edge(51), Grid::EdgeNodes({{23, 26}}));
    EXPECT_EQ(grid.edge(52), Grid::EdgeNodes({{24, 25}}));
    EXPECT_EQ(grid.edge(53), Grid::EdgeNodes({{25, 26}}));

    // Cell Node indices
    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_cells(); ++i) {
        const auto node_indices = grid.node_indices_of((Grid::CellIndex) i);
        const Hexa h(
            grid.node(node_indices[0]), grid.node(node_indices[1]), grid.node(node_indices[2]),
            grid.node(node_indices[3]),
            grid.node(node_indices[4]), grid.node(node_indices[5]), grid.node(node_indices[6]),
            grid.node(node_indices[7]));
        EXPECT_MATRIX_NEAR(h.nodes(), grid.cell_at((Grid::CellIndex) i).nodes(), 1e-15);
    }

    // Face queries
    EXPECT_EQ(grid.number_of_faces(), (unsigned) 36);
    // First slice (2D grid)
    EXPECT_EQ(grid.face(0), Grid::FaceNodes({{0, 1, 4, 3}}));
    EXPECT_EQ(grid.face(1), Grid::FaceNodes({{1, 2, 5, 4}}));
    EXPECT_EQ(grid.face(2), Grid::FaceNodes({{3, 4, 7, 6}}));
    EXPECT_EQ(grid.face(3), Grid::FaceNodes({{4, 5, 8, 7}}));
    // Between the slice 1 and slice 2
    // - xz axis
    EXPECT_EQ(grid.face(4), Grid::FaceNodes({{0, 1, 10,  9}}));
    EXPECT_EQ(grid.face(5), Grid::FaceNodes({{1, 2, 11, 10}}));
    EXPECT_EQ(grid.face(6), Grid::FaceNodes({{3, 4, 13, 12}}));
    EXPECT_EQ(grid.face(7), Grid::FaceNodes({{4, 5, 14, 13}}));
    EXPECT_EQ(grid.face(8), Grid::FaceNodes({{6, 7, 16, 15}}));
    EXPECT_EQ(grid.face(9), Grid::FaceNodes({{7, 8, 17, 16}}));
    // - yz axis
    EXPECT_EQ(grid.face(10), Grid::FaceNodes({{0, 3, 12,  9}}));
    EXPECT_EQ(grid.face(11), Grid::FaceNodes({{1, 4, 13, 10}}));
    EXPECT_EQ(grid.face(12), Grid::FaceNodes({{2, 5, 14, 11}}));
    EXPECT_EQ(grid.face(13), Grid::FaceNodes({{3, 6, 15, 12}}));
    EXPECT_EQ(grid.face(14), Grid::FaceNodes({{4, 7, 16, 13}}));
    EXPECT_EQ(grid.face(15), Grid::FaceNodes({{5, 8, 17, 14}}));

    // Second slice (2D grid)
    int s = 9; // slice_corner_node_index
    EXPECT_EQ(grid.face(16), Grid::FaceNodes({{s + 0, s + 1, s + 4, s + 3}}));
    EXPECT_EQ(grid.face(17), Grid::FaceNodes({{s + 1, s + 2, s + 5, s + 4}}));
    EXPECT_EQ(grid.face(18), Grid::FaceNodes({{s + 3, s + 4, s + 7, s + 6}}));
    EXPECT_EQ(grid.face(19), Grid::FaceNodes({{s + 4, s + 5, s + 8, s + 7}}));
    // Between the slice 2 and slice 3
    // - xz axis
    EXPECT_EQ(grid.face(20), Grid::FaceNodes({{s + 0, s + 1, s + 10, s + 9}}));
    EXPECT_EQ(grid.face(21), Grid::FaceNodes({{s + 1, s + 2, s + 11, s + 10}}));
    EXPECT_EQ(grid.face(22), Grid::FaceNodes({{s + 3, s + 4, s + 13, s + 12}}));
    EXPECT_EQ(grid.face(23), Grid::FaceNodes({{s + 4, s + 5, s + 14, s + 13}}));
    EXPECT_EQ(grid.face(24), Grid::FaceNodes({{s + 6, s + 7, s + 16, s + 15}}));
    EXPECT_EQ(grid.face(25), Grid::FaceNodes({{s + 7, s + 8, s + 17, s + 16}}));
    // - yz axis
    EXPECT_EQ(grid.face(26), Grid::FaceNodes({{s + 0, s + 3, s + 12, s + 9}}));
    EXPECT_EQ(grid.face(27), Grid::FaceNodes({{s + 1, s + 4, s + 13, s + 10}}));
    EXPECT_EQ(grid.face(28), Grid::FaceNodes({{s + 2, s + 5, s + 14, s + 11}}));
    EXPECT_EQ(grid.face(29), Grid::FaceNodes({{s + 3, s + 6, s + 15, s + 12}}));
    EXPECT_EQ(grid.face(30), Grid::FaceNodes({{s + 4, s + 7, s + 16, s + 13}}));
    EXPECT_EQ(grid.face(31), Grid::FaceNodes({{s + 5, s + 8, s + 17, s + 14}}));

    // Third slice (2D grid)
    s = 18; // slice_corner_node_index
    EXPECT_EQ(grid.face(32), Grid::FaceNodes({{s + 0, s + 1, s + 4, s + 3}}));
    EXPECT_EQ(grid.face(33), Grid::FaceNodes({{s + 1, s + 2, s + 5, s + 4}}));
    EXPECT_EQ(grid.face(34), Grid::FaceNodes({{s + 3, s + 4, s + 7, s + 6}}));
    EXPECT_EQ(grid.face(35), Grid::FaceNodes({{s + 4, s + 5, s + 8, s + 7}}));

    // Cell queries by position
    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_cells(); ++i) {
        const auto cell = grid.cell_at(i);
        for (const auto & gauss_node : cell.gauss_nodes()) {
            const auto p = cell.world_coordinates(gauss_node.position);
            EXPECT_EQ(grid.cell_index_containing(p, false), i);
        }
    }
}

TEST(Topology_Grid_3D, BenchMark)
{
    using namespace caribou::topology;
    using Grid = Grid<3>;

    using WorldCoordinates = Grid::WorldCoordinates;
    using Subdivisions = Grid::Subdivisions;
    using Dimensions = Grid::Dimensions;
    BEGIN_CLOCK;

    Grid grid(WorldCoordinates{0, 0, 0}, Subdivisions{10, 10, 10}, Dimensions{100, 100, 100});

    std::vector<WorldCoordinates> positions(grid.number_of_nodes());
    std::vector<Grid::ElementNodes> hexahedrons(grid.number_of_cells());

    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_nodes(); ++i) {
        positions[i] = grid.node(i);
    }

    for (UNSIGNED_INTEGER_TYPE i = 0; i < grid.number_of_cells(); ++i) {
        hexahedrons[i] = grid.node_indices_of(i);
    }

    // Memory-bound test
    FLOATING_POINT_TYPE volume_1 = 0;
    TICK;
    for (UNSIGNED_INTEGER_TYPE i = 0; i < hexahedrons.size(); ++i) {
        const auto & node_indices = hexahedrons[i];
        std::array<WorldCoordinates, 8> nodes;
        for (UNSIGNED_INTEGER_TYPE j = 0; j < 8; ++j) {
            nodes[j] = positions[node_indices[j]];
        }

        for (UNSIGNED_INTEGER_TYPE n = 0; n < 8; ++n) {
            for (UNSIGNED_INTEGER_TYPE m = 0; m < 8; ++m) {
                volume_1 += (nodes[n] - nodes[m]).norm();
            }
        }
    }
    std::cout <<"Memory-bound grid 3d: " << TOCK / 1000. / 1000. << " [ms]" << std::endl;


    // CPU-bound test
    FLOATING_POINT_TYPE volume_2 = 0;
    const auto ncells = grid.number_of_cells();
    TICK;
    for (UNSIGNED_INTEGER_TYPE i = 0; i < ncells; ++i) {
        Grid::ElementNodes node_indices = grid.node_indices_of(i);
        std::array<WorldCoordinates, 8> nodes;
        for (UNSIGNED_INTEGER_TYPE j = 0; j < 8; ++j) {
            nodes[j] = grid.node(node_indices[j]);
        }

        for (UNSIGNED_INTEGER_TYPE n = 0; n < 8; ++n) {
            for (UNSIGNED_INTEGER_TYPE m = 0; m < 8; ++m) {
                volume_2 += (nodes[n] - nodes[m]).norm();
            }
        }
    }
    std::cout <<"CPU-bound grid 3d: " << TOCK / 1000. / 1000. << " [ms]" << std::endl;

    EXPECT_EQ(volume_1, volume_2);
}

#endif //CARIBOU_TOPOLOGY_TEST_GRID_3D_H
