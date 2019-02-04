#include <gtest/gtest.h>

#include <Caribou/Algebra/Matrix.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>


static constexpr double youngModulus = 100000;
static constexpr double poissonRatio = 0;

static constexpr double l = youngModulus*poissonRatio / ((1 + poissonRatio)*(1 - 2*poissonRatio));
static constexpr double m = youngModulus / (2 * (1 + poissonRatio));
static constexpr double a = l + 2*m;
static const caribou::algebra::Matrix<6,6,double> C ({
    {a,l,l,0,0,0},
    {l,a,l,0,0,0},
    {l,l,a,0,0,0},
    {0,0,0,m,0,0},
    {0,0,0,0,m,0},
    {0,0,0,0,0,m}
});

TEST(Mechanics, Strain) {
    using namespace caribou::geometry;
    using namespace caribou::algebra;
    using namespace caribou::mechanics;

    using LocalCoordinates = Vector<3, FLOATING_POINT_TYPE>;

    using Mat33 = Matrix<3,3>;
    const auto I = Mat33::Identity();

    Hexahedron<interpolation::Hexahedron8> initial_hexahedron;
    {
        const Mat33 F = elasticity::strain::F(initial_hexahedron, initial_hexahedron, LocalCoordinates {1/sqrt(3), 1/sqrt(3), 1/sqrt(3)});
        ASSERT_EQ(Mat33::Identity(), F);
    }

    Matrix<8,3> displacements ({
        {-0.5, -1.2, -0.8},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {-0.6,  0.7,  0.7}
    });

    Vector<24> U;
    for (size_t i = 0; i < 8; ++i) {
        U[i * 3 + 0] = displacements(i,0);
        U[i * 3 + 1] = displacements(i,1);
        U[i * 3 + 2] = displacements(i,2);
    }

    Hexahedron<interpolation::Hexahedron8> deformed_hexahedron = initial_hexahedron;
    for (size_t i = 0; i < 8; ++i)
        deformed_hexahedron.node(i) += displacements.row(i).T();

    // Compute the elastic force with the stiffness matrix
    Matrix<24,24,double> K = initial_hexahedron.gauss_quadrature([](const auto & hexa, const auto & local_coordinates) {
        const auto B = elasticity::strain::B(hexa, local_coordinates);
        return B.T() * C * B;
    });
    Vector<24> F1 = K*U;

    // Compute the elastic force by manually integrating the stress and strain tensors
    Vector<24> F2; F2.fill(0);
    for (const auto & gauss_node : interpolation::Hexahedron8::gauss_nodes) {
        const auto J = initial_hexahedron.jacobian(gauss_node);
        const auto Jinv = J^-1;
        const auto detJ = J.determinant();

        Matrix<8,3> dN_dx = (Jinv.T() * interpolation::Hexahedron8::dN(gauss_node).T()).T();

        const auto F = elasticity::strain::F(dN_dx, displacements);
        const auto E2 = (F.T() + F) - 2*I;
        const auto S = m*E2 + 0.5*(l*tr(E2)*I);

        for (size_t i = 0; i < 8; ++i) {
            const auto dx = dN_dx.row(i).T();
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
    const auto Jinv = J^-1;
    const auto detJ = J.determinant();

    for (const auto & gauss_node : interpolation::Hexahedron8::gauss_nodes) {

        Matrix<8,3> dN_dx = (Jinv.T() * interpolation::Hexahedron8::dN(gauss_node).T()).T();

        const auto F = elasticity::strain::F(dN_dx, displacements);
        const auto E2 = (F.T() + F) - 2*I;
        const auto S = m*E2 + 0.5*(l*tr(E2)*I);

        for (size_t i = 0; i < 8; ++i) {
            const auto dx = dN_dx.row(i).T();
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
