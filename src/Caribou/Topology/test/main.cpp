#include <gtest/gtest.h>
#include <Caribou/Topology/Grid.h>

TEST(Topology, Grid2D) {
    using namespace caribou::topology;

    Grid2D grid(
            {0, 0} /* anchor */,
            {4, 4} /* subdivision */,
            {100, 100} /* size */
    );

    ASSERT_THROW(grid.get({4,4}), std::out_of_range);
    ASSERT_FLOAT_EQ(100/4., grid.get({0,0}).size()[0]);
    ASSERT_EQ((Grid2D::Index)0, grid.get({0,0}).index());
    ASSERT_EQ((Grid2D::Index)15, grid.get({3,3}).index());
    ASSERT_EQ((Grid2D::Index)24, grid.nodes({3,3})[3]);
    ASSERT_EQ(Grid2D::VecInt ({2,1,2}), grid.grid_coordinates(grid.cell_index({2,1,2})));


//    auto H = grid.cell_size();
//    const auto hx = H[0];
//    ASSERT_FLOAT_EQ(100/9., hx);
//
//    Grid & cell = grid.query(4, 4, 4);
//    ASSERT_TRUE(cell.is_a_leaf());
//    cell.subdivide();
//    ASSERT_FALSE(cell.is_a_leaf());
//
//    H = cell.cell_size();
//    ASSERT_FLOAT_EQ(hx/2., H[0]);
//
//    VecFloat middle = cell.position({0, 0, 0});
//    ASSERT_EQ(VecFloat(50, 50, 50), middle);
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
