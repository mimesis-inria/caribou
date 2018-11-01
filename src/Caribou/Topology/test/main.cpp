#include <gtest/gtest.h>
#include <Caribou/Topology/Grid.h>

TEST(Geometry, Topology) {
    using namespace caribou::topology;

    using VecFloat = Grid::VecFloat;

    Grid grid({0,0,0} /* anchor */, {9, 9, 9} /* subdivision */, {100, 100, 100} /* size */);

    auto H = grid.cell_size();
    const auto hx = H[0];
    ASSERT_FLOAT_EQ(100/9., hx);

    Grid & cell = grid.query(4, 4, 4);
    ASSERT_TRUE(cell.is_a_leaf());
    cell.subdivide();
    ASSERT_FALSE(cell.is_a_leaf());

    H = cell.cell_size();
    ASSERT_FLOAT_EQ(hx/2., H[0]);

    VecFloat middle = cell.position({0, 0, 0});
    ASSERT_EQ(VecFloat(50, 50, 50), middle);
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
