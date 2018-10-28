#include <gtest/gtest.h>
#include <Caribou/Algebra/Vector.h>

TEST(Geometry, Vector) {
    using namespace caribou::algebra;

    // Constructor and assigments tests
    Vector<3> v1 {1, 2, 3};
    Vector<3> v2 ({1,2,3});
    Vector<3> v3 (1, 2, 3);
    Vector<3> v4 (0,0,0);
    Vector<3> v5(true); // Automatically fills the vector with 0 values


    ASSERT_EQ(v1, v2);
    ASSERT_EQ(v2, v3);
    ASSERT_EQ(v4, v5);

    Vector<3, char> v6 ( (float) 0, (float) 0, (float) 0);
    ASSERT_EQ(typeid(v6[2]), typeid(char));
    ASSERT_EQ(sizeof(v6), 3*sizeof(char));

    Vector<3, float> v7 = {(char) 1, 2, 3};
    ASSERT_EQ(v1, v7);
    ASSERT_EQ(typeid(v7[2]), typeid(float));

    auto scalar = v2 * v3;
    ASSERT_EQ(scalar, v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2]);

    auto mult = v2*3;
    ASSERT_EQ(mult, Vector<3>(v2[0]*3, v2[1]*3, v2[2]*3));

    auto sum = v2 + v3;
    ASSERT_EQ(sum, Vector<3>(v2[0]+v3[0], v2[1]+v3[1], v2[2]+v3[2]));

    // Cross product
    Vector<3> v8 {10, 0, 0};
    Vector<3> v9 {0, 10, 0};
    Vector<3> v10 {0, 0, 100};
    ASSERT_EQ(v8^v9, v10);

    // Unitary vector
    Vector<3> v11 {0, 0, 1};
    ASSERT_EQ(v10.unit(), v11);
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
