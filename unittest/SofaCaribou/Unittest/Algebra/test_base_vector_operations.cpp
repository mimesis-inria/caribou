#include <gtest/gtest.h>

#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <SofaBaseLinearSolver/FullVector.h>
#include <SofaCaribou/Algebra/BaseVectorOperations.h>
DISABLE_ALL_WARNINGS_END

#include <Eigen/Dense>

TEST(Algebra, SofaFullDFullDDotProduct) {
    const auto n = 100;
    const Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
    const Eigen::VectorXd v2 = Eigen::VectorXd::Random(n);

    sofa::component::linearsolver::FullVector<double> sofa_v1 (n);
    sofa::component::linearsolver::FullVector<double> sofa_v2 (n);

    for (sofa::Index i = 0; i < n; ++i) {
        sofa_v1[i] = v1[static_cast<Eigen::Index>(i)];
        sofa_v2[i] = v2[static_cast<Eigen::Index>(i)];
    }

    EXPECT_NEAR(SofaCaribou::Algebra::dot(&sofa_v1, &sofa_v2), v1.dot(v2), 1e-10);
}

TEST(Algebra, SofaFullDFullFDotProduct) {
    const auto n = 100;
    const Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
    const Eigen::VectorXf v2 = Eigen::VectorXf::Random(n);

    sofa::component::linearsolver::FullVector<double> sofa_v1 (n);
    sofa::component::linearsolver::FullVector<float> sofa_v2 (n);

    for (sofa::Index i = 0; i < n; ++i) {
        sofa_v1[i] = v1[static_cast<Eigen::Index>(i)];
        sofa_v2[i] = v2[static_cast<Eigen::Index>(i)];
    }

    EXPECT_NEAR(SofaCaribou::Algebra::dot(&sofa_v1, &sofa_v2), v1.dot(v2.cast<double>()), 1e-10);
}

TEST(Algebra, SofaFullFFullDDotProduct) {
    const auto n = 100;
    const Eigen::VectorXf v1 = Eigen::VectorXf::Random(n);
    const Eigen::VectorXd v2 = Eigen::VectorXd::Random(n);

    sofa::component::linearsolver::FullVector<float> sofa_v1 (n);
    sofa::component::linearsolver::FullVector<double> sofa_v2 (n);

    for (sofa::Index i = 0; i < n; ++i) {
        sofa_v1[i] = v1[static_cast<Eigen::Index>(i)];
        sofa_v2[i] = v2[static_cast<Eigen::Index>(i)];
    }

    EXPECT_NEAR(SofaCaribou::Algebra::dot(&sofa_v1, &sofa_v2), v1.cast<double>().dot(v2), 1e-10);
}