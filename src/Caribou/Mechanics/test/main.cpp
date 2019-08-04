#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

template<int nRows, int nColumns, int Options=0>
using Matrix = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, nColumns, Options>;

template<int nRows, int Options=0>
using Vector = Eigen::Matrix<FLOATING_POINT_TYPE, nRows, 1, Options>;

template<int nRows>
using MapVector = Eigen::Map<const Vector<nRows, Eigen::ColMajor>>;

static constexpr double youngModulus = 100000;
static constexpr double poissonRatio = 0;

static constexpr double l = youngModulus*poissonRatio / ((1 + poissonRatio)*(1 - 2*poissonRatio));
static constexpr double m = youngModulus / (2 * (1 + poissonRatio));
static constexpr double a = l + 2*m;


TEST(Mechanics, Strain) {
    using namespace caribou::geometry;
    using namespace caribou::mechanics;

    using LocalCoordinates = Vector <3>;

    using Mat33 = Matrix<3,3>;
    const auto I = Mat33::Identity();

    Matrix<6,6> C;
    C <<
        a,l,l,0,0,0,
        l,a,l,0,0,0,
        l,l,a,0,0,0,
        0,0,0,m,0,0,
        0,0,0,0,m,0,
        0,0,0,0,0,m;

    Hexahedron<interpolation::Hexahedron8> initial_hexahedron;
    {
        const Mat33 F = elasticity::strain::F(initial_hexahedron, initial_hexahedron, LocalCoordinates {1/sqrt(3), 1/sqrt(3), 1/sqrt(3)});
        ASSERT_EQ(Mat33::Identity(), F);
    }

    Matrix<8,3> displacements;
    displacements <<
        -0.5, -1.2, -0.8,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        -0.6,  0.7,  0.7;

    Vector<24> U;
    for (size_t i = 0; i < 8; ++i) {
        U[i * 3 + 0] = displacements(i,0);
        U[i * 3 + 1] = displacements(i,1);
        U[i * 3 + 2] = displacements(i,2);
    }

    Hexahedron<interpolation::Hexahedron8> deformed_hexahedron = initial_hexahedron;
    for (size_t i = 0; i < 8; ++i)
        deformed_hexahedron.node(i) += displacements.row(i).transpose();

    // Compute the elastic force with the stiffness matrix
    auto K = initial_hexahedron.gauss_quadrature<Matrix<24,24>>([&C](const auto & hexa, const auto & local_coordinates) {
        const auto B = elasticity::strain::B(hexa, local_coordinates);
        const Matrix<24, 24> K = B.transpose() * C * B;
        return K;
    });
    Vector<24> F1 = K*U;

    // Compute the elastic force by manually integrating the stress and strain tensors
    Vector<24> F2; F2.fill(0);
    for (std::size_t gauss_node_id = 0; gauss_node_id < interpolation::Hexahedron8::number_of_gauss_nodes; ++gauss_node_id) {
        const auto &gauss_node   = MapVector<3>(interpolation::Hexahedron8::gauss_nodes[gauss_node_id]);

        const auto J = initial_hexahedron.jacobian(gauss_node);
        const auto Jinv = J.inverse();
        const auto detJ = J.determinant();

        Matrix<8,3> dN_dx = (Jinv.transpose() * interpolation::Hexahedron8::dL(gauss_node).transpose()).transpose();

        const auto F = elasticity::strain::F(dN_dx, displacements);
        const auto E2 = (F.transpose() + F) - 2*I;
        const auto S = m*E2 + 0.5*(l*E2.trace()*I);

        for (size_t i = 0; i < 8; ++i) {
            const auto dx = dN_dx.row(i).transpose();
            const auto f = S*dx*detJ;
            F2[i*3 + 0] += f[0];
            F2[i*3 + 1] += f[1];
            F2[i*3 + 2] += f[2];
        }
    }

    for (size_t i = 0; i < 24; ++i)
        ASSERT_FLOAT_EQ(F1[i], F2[i]);

    // Same thing, but taking into account the rectangularity of the hexa
    F2.fill(0);
    const auto J = Matrix<3,3>::Identity();
    const auto Jinv = J.inverse();
    const auto detJ = J.determinant();

    for (std::size_t gauss_node_id = 0; gauss_node_id < interpolation::Hexahedron8::number_of_gauss_nodes; ++gauss_node_id) {
        const auto &gauss_node   = MapVector<3>(interpolation::Hexahedron8::gauss_nodes[gauss_node_id]);

        Matrix<8,3> dN_dx = (Jinv.transpose() * interpolation::Hexahedron8::dL(gauss_node).transpose()).transpose();

        const auto F = elasticity::strain::F(dN_dx, displacements);
        const auto E2 = (F.transpose() + F) - 2*I;
        const auto S = m*E2 + 0.5*(l*E2.trace()*I);

        for (size_t i = 0; i < 8; ++i) {
            const auto dx = dN_dx.row(i).transpose();
            const auto f = S*dx*detJ;
            F2[i*3 + 0] += f[0];
            F2[i*3 + 1] += f[1];
            F2[i*3 + 2] += f[2];
        }
    }

    for (size_t i = 0; i < 24; ++i)
        ASSERT_FLOAT_EQ(F1[i], F2[i]);

}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
