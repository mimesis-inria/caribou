#include <gtest/gtest.h>

#include <SofaBaseLinearSolver/FullVector.h>
#include <SofaCaribou/Algebra/BaseVectorOperations.h>
#include <Caribou/config.h>
#include <Eigen/Dense>

TEST(Algebra, SofaFullDFullDDotProduct) {
    const auto n = 100;
    const Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
    const Eigen::VectorXd v2 = Eigen::VectorXd::Random(n);

    sofa::component::linearsolver::FullVector<double> sofa_v1 (n);
    sofa::component::linearsolver::FullVector<double> sofa_v2 (n);

    for (std::size_t i = 0; i < n; ++i) {
        sofa_v1[i] = v1[i];
        sofa_v2[i] = v2[i];
    }

    EXPECT_NEAR(SofaCaribou::Algebra::dot(&sofa_v1, &sofa_v2), v1.dot(v2), 1e-10);
}

TEST(Algebra, SofaFullDFullFDotProduct) {
    const auto n = 100;
    const Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
    const Eigen::VectorXd v2 = Eigen::VectorXf::Random(n).cast<double>();

    sofa::component::linearsolver::FullVector<double> sofa_v1 (n);
    sofa::component::linearsolver::FullVector<float> sofa_v2 (n);

    for (std::size_t i = 0; i < n; ++i) {
        sofa_v1[i] = v1[i];
        sofa_v2[i] = v2[i];
    }

    EXPECT_NEAR(SofaCaribou::Algebra::dot(&sofa_v1, &sofa_v2), v1.dot(v2), 1e-10);
}

TEST(Algebra, SofaFullFFullDDotProduct) {
    const auto n = 100;
    const Eigen::VectorXd v1 = Eigen::VectorXf::Random(n).cast<double>();
    const Eigen::VectorXd v2 = Eigen::VectorXd::Random(n);

    sofa::component::linearsolver::FullVector<float> sofa_v1 (n);
    sofa::component::linearsolver::FullVector<double> sofa_v2 (n);

    for (std::size_t i = 0; i < n; ++i) {
        sofa_v1[i] = v1[i];
        sofa_v2[i] = v2[i];
    }

    EXPECT_NEAR(SofaCaribou::Algebra::dot(&sofa_v1, &sofa_v2), v1.dot(v2), 1e-10);
}