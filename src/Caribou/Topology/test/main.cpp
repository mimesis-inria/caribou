#include <gtest/gtest.h>
#include <Caribou/Topology/Engine/Grid/Grid.h>

TEST(Topology, Grid2D) {
    using namespace caribou::topology::engine;

    Grid2D grid(
            {0.75, 0.5} /* anchor */,
            {4, 4} /* subdivision */,
            {100, 100} /* size */
    );

    ASSERT_THROW(grid.get({4,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,4}), std::out_of_range);

    ASSERT_FLOAT_EQ(100/4., grid.get({0,0}).size()[0]);

    ASSERT_EQ((Grid2D::Index)0, grid.get({0,0}).index());
    ASSERT_EQ((Grid2D::Index)15, grid.get({3,3}).index());
    ASSERT_EQ((Grid2D::Index)24, grid.nodes({3,3})[3]);
    ASSERT_EQ(Grid2D::VecInt ({2,1,2}), grid.grid_coordinates(grid.cell_index({2,1,2})));

    using Cell = Grid2D::CellType;

    const Cell & cell = grid.get({1,2});
    const Grid2D::Index node_id = cell.nodes()[0];
    ASSERT_EQ((Grid2D::Index)11, node_id);
    ASSERT_EQ(Grid2D::VecFloat({25.75, 50.5}), grid.position(node_id));
}

TEST(Topology, Grid3D) {
    using namespace caribou::topology::engine;
    Grid3D grid(
            {0, 0, 0} /* anchor */,
            {2, 2, 2} /* subdivision */,
            {100, 100, 100} /* size */
    );

    ASSERT_THROW(grid.get({2,0,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,2,0}), std::out_of_range);
    ASSERT_THROW(grid.get({0,0,2}), std::out_of_range);

    ASSERT_FLOAT_EQ(100/2., grid.get({1,0,1}).size()[0]);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
