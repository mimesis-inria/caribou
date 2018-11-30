#include <gtest/gtest.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Algebra/Matrix.h>

TEST(Algebra, Vector) {
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
    Matrix<1, 3, char> v8 = {(char) 1, 2, 3};
    ASSERT_EQ(v1, v7);
    ASSERT_EQ(v1, v8);
    ASSERT_EQ(typeid(v7[2]), typeid(float));
    ASSERT_EQ(typeid(v8[2]), typeid(char));

    // Inner product (dot or scalar product)
    auto scalar = v2.inner_product(v3);
    ASSERT_EQ(scalar, v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2]);

    scalar = v2.inner_product(v8);
    ASSERT_EQ(scalar, v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2]);

    scalar = v8 * v2;
    ASSERT_EQ(scalar, v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2]);

    // Scalar multiplication
    auto mult = v2*3;
    ASSERT_EQ(mult, Vector<3>(v2[0]*3, v2[1]*3, v2[2]*3));

    mult = 3*v2;
    ASSERT_EQ(mult, Vector<3>(v2[0]*3, v2[1]*3, v2[2]*3));

    v2 *= 3;
    ASSERT_EQ(v2, mult);

    v2 /= 3;
    ASSERT_EQ(v2, v1);

    auto sum = v2 + v3;
    ASSERT_EQ(sum, Vector<3>(v2[0]+v3[0], v2[1]+v3[1], v2[2]+v3[2]));

    // Cross product
    Vector<3> v9 {10, 0, 0};
    Vector<3> v10 {0, 10, 0};
    Vector<3> v11 {0, 0, 100};
    ASSERT_EQ(v9.cross(v10), v11);
    ASSERT_EQ(v9^v10, v11);

    // Direct sum
    v11 =  Vector<3>(10, 10, 0);
    ASSERT_EQ(v9.direct_sum(v10), v11);

    // Direct sub
    v11 =  Vector<3>(10, -10, 0);
    ASSERT_EQ(v9.direct_sub(v10), v11);

    // Direct mult
    v11 =  Vector<3>(0, 0, 0);
    ASSERT_EQ(v9.direct_mult(v10), v11);

    // Unitary vector
    Vector<3> v12 {0, 1, 0};
    ASSERT_EQ(v10.unit(), v12);
}

TEST(Algebra, Matrix) {
    using namespace caribou::algebra;

    // Constructor
    Matrix<4,3> m1 (
        {
            {1,2,3},  // Row 0
            {4,5,6},  // Row 1
            {7,8,9},  // Row 2
            {10,11,12}  // Row 3
        }
    );
    ASSERT_EQ(8, m1(2,1));

    Matrix<3,4> m1_t = m1.transposed();
    ASSERT_EQ(m1(1,2), m1_t(2,1));

    // Identity
    Matrix<4,4> I = {{{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}}};
    auto I_ = Matrix<4,4>::Identity();
    ASSERT_EQ(I, I_);

    // Determinant
    Matrix<2,2> m2_2 ({{5,2}, {9,4}});
    Matrix<3,3> m2_3 (
            {
                    { 3,0, 2},  // Row 0
                    {-2,0,-2},  // Row 1
                    { 0,1, 1}  // Row 2
            }
    );

    ASSERT_EQ(2, m2_2.determinant());
    ASSERT_EQ(2, m2_3.determinant());

    // Inverse
    Matrix<2,2> m2_2_inverted (
            {
                    { 2, -1},  // Row 0
                    {-9/2., 5/2.}  // Row 1
            }
    );
    Matrix<3,3> m2_3_inverted (
            {
                    { 1, 1, 0},  // Row 0
                    {1, 3/2., 1},  // Row 1
                    {-1, -3/2., 0}  // Row 2
            }
    );

    ASSERT_EQ(m2_2_inverted, m2_2.inverted());
    ASSERT_EQ(m2_2_inverted, m2_2^-1);
    ASSERT_EQ(m2_3_inverted, m2_3.inverted());
    ASSERT_EQ(m2_3_inverted, m2_3^-1);

    // Operation
    Matrix<3,3> m3 (
            {
                    {1,2,3},  // Row 0
                    {4,5,6},  // Row 1
                    {7,8,9},  // Row 2
            }
    );

    Matrix<3,1> v1 ({
        {10},
        {11},
        {12}
    });

    Vector<3> v2 {10, 11, 12};

    Matrix<3,1> res ({
        {10*1 + 11*2 + 12*3},
        {10*4 + 11*5 + 12*6},
        {10*7 + 11*8 + 12*9}
    });

    ASSERT_EQ(res, m3*v1);
    ASSERT_EQ(res, m3*v2);

    Matrix<3,3> res2 = {{
        {2,4,6},  // Row 0
        {8,10,12},  // Row 1
        {14,16,18},  // Row 2
    }};

    Matrix<3,3> m4 = {{
        {1,2,3},  // Row 0
        {4,5,6},  // Row 1
        {7,8,9},  // Row 2
    }};

    ASSERT_EQ(res2, m4+m4);
    m4 += m4;
    ASSERT_EQ(res2, m4);

    Vector <3> v3 {1, 2, 3};
    Matrix<3,3> res3 {{
        {1, 2, 3},
        {2, 4, 6},
        {3, 6, 9}
    }};
    ASSERT_EQ(res3, v3 * v3.T());


}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
