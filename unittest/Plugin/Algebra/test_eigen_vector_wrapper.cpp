#include <gtest/gtest.h>
#include <SofaCaribou/Algebra/EigenVector.h>

#include <Caribou/config.h>
#include <Eigen/Dense>

template<int nRows, int nColumns>
using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns>;

template <typename Derived>
void test_core(SofaCaribou::Algebra::EigenVector<Derived> & mm) {
    using Index = typename SofaCaribou::Algebra::EigenVector<Derived>::Index;

    // Set-get
    mm.set(30, 200);
    EXPECT_EQ(mm.element(30), 200);

    // Add
    mm.add(30, -100);
    EXPECT_EQ(mm.element(30), 100);

    // Resize
    mm.resize(50);
    EXPECT_EQ(mm.size(), 50);
    int nb_not_equal = 0;
    for (Index i=0;i<mm.size();++i)
            if (mm.element(i) != 0) nb_not_equal++;

    EXPECT_EQ(nb_not_equal, 0) << "There are " << nb_not_equal << " values that are not equal to zero (and they should).";

}

TEST(Algebra, ColVectorCreation) {
    using Vector = Eigen::Matrix<float, Eigen::Dynamic, 1>;
    using EigenVector = SofaCaribou::Algebra::EigenVector<Vector>;
    EigenVector mm(100);
    EXPECT_EQ(mm.size(), 100);

    // Testing the core functionalities
    test_core(mm);
}

TEST(Algebra, ColVectorByCopy) {
    using Vector = Eigen::Matrix<float, Eigen::Dynamic, 1>;
    Vector m(100);

    using EigenVector = SofaCaribou::Algebra::EigenVector<Vector>;
    EigenVector mm(m);
    EXPECT_EQ(mm.size(), 100);

    // Testing the core functionalities
    test_core(mm);

    // Set - should not change the values of the initial Eigen matrix since the wrapper made a copy
    m[30] = 100;
    mm.set(30, 200);
    EXPECT_EQ(mm.element(30), 200);
    EXPECT_EQ(m[30], 100);
}

TEST(Algebra, ColVectorByReference) {
    using Vector = Eigen::Matrix<double, 1, Eigen::Dynamic>;
    using EigenVector = SofaCaribou::Algebra::EigenVector<Vector &>;
    Vector m(100);

    EigenVector mm(m);

    EXPECT_EQ(mm.size(), 100);

    // Testing the core functionalities
    test_core(mm);

    // Set - should change the values of the initial Eigen matrix since the wrapper has a reference.
    m[30] = 100;
    mm.set(30, 200);
    EXPECT_EQ(mm.element(30), 200);
    EXPECT_EQ(m[30], 200);
}
