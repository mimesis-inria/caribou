#include <gtest/gtest.h>
#include <Caribou/Topology/Engine/Grid/Grid.h>

TEST(Topology, Grid2D) {
    using namespace caribou::topology::engine;
    using Cell = Grid2D::CellType;

    Grid2D grid(
            {0.75, 0.5} /* anchor */,
            {4, 4} /* subdivision */,
            {100, 100} /* size */
    );

    // Grid's cell accessor test
    ASSERT_THROW(grid.get({4,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,4}), std::out_of_range);

    // Cell size test
    ASSERT_FLOAT_EQ(100/4., grid.get({0,0}).size()[0]);

    // Node numbering test
    ASSERT_EQ((Grid2D::Index)0, grid.get({0,0}).index());
    ASSERT_EQ((Grid2D::Index)15, grid.get({3,3}).index());
    ASSERT_EQ((Grid2D::Index)24, grid.nodes({3,3})[2]);
    ASSERT_EQ(Grid2D::VecInt ({2,1,2}), grid.grid_coordinates(grid.cell_index({2,1,2})));

    // Positioning test
    Cell & cell = grid.get({1,2});
    Grid2D::Index node_id = cell.nodes()[0];
    ASSERT_EQ((Grid2D::Index)11, node_id);
    ASSERT_EQ(Grid2D::VecFloat({25.75, 50.5}), grid.position(node_id));

    // Subdivision test
    Grid2D & sub_grid = *cell.subdivide({2,2});
    Cell & sub_cell = sub_grid.get({1, 0});
    node_id = sub_cell.nodes()[1];
    ASSERT_EQ((Grid2D::Index)2, node_id);
    ASSERT_EQ(Grid2D::VecFloat({50.75, 50.5}), sub_grid.position(node_id));
}

TEST(Topology, Grid3D) {
    using namespace caribou::topology::engine;
    using Cell = Grid3D::CellType;

    Grid3D grid(
            {0.25, 0.5, 0.75} /* anchor */,
            {2, 2, 2} /* subdivision */,
            {100, 100, 100} /* size */
    );

    // Grid's cell accessor test
    ASSERT_THROW(grid.get({2,0,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,2,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,0,2}), std::out_of_range);

    // Cell size test
    ASSERT_FLOAT_EQ(100/2., grid.get({1,0,1}).size()[0]);

    // Node numbering test
    ASSERT_EQ((Grid3D::Index)0, grid.get({0,0,0}).index());
    ASSERT_EQ((Grid3D::Index)3, grid.get({1,1,0}).index());
    ASSERT_EQ((Grid3D::Index)6, grid.get({0,1,1}).index());

    // Positioning test
    Cell & cell = grid.get({1,0,1});
    Grid3D::Index node_id = cell.nodes()[0];
    ASSERT_EQ((Grid3D::Index)10, node_id);
    ASSERT_EQ(Grid3D::VecFloat({50.25, 0.5, 50.75}), grid.position(node_id));

    // Subdivision test
    Grid3D & sub_grid = *cell.subdivide({2,2,2});
    Cell & sub_cell = sub_grid.get({0,1,1});
    node_id = sub_cell.nodes()[7];
    ASSERT_EQ((Grid3D::Index)24, node_id);
    ASSERT_EQ(Grid3D::VecFloat({50.25, 50.5, 100.75}), sub_grid.position(node_id));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
