#pragma once

#include <Caribou/Algebra/Tensor.h>
#include <SofaCaribou/Material/HyperelasticMaterial.h>

namespace SofaCaribou::material {

template<class DataTypes>
class MooneyRivlinMaterial : public HyperelasticMaterial<DataTypes> {
    static constexpr auto Dimension = DataTypes::spatial_dimensions;
    using Coord = typename DataTypes::Coord;
    using Real  = typename Coord::value_type;
public:
    SOFA_CLASS(SOFA_TEMPLATE(MooneyRivlinMaterial, DataTypes), SOFA_TEMPLATE(HyperelasticMaterial, DataTypes));

    MooneyRivlinMaterial()
        : d_c01(initData(&d_c01,
            Real(1000), "c01",
            "Young's modulus of the material",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_c10(initData(&d_c10,
            Real(0.3),  "c10",
           "Poisson's ratio of the material",
           true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_k(initData(&d_k,
            Real(0.3),  "k",
           "Poisson's ratio of the material",
           true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
    {
    }

    /**
     * This is called just before the material is updated on every points (usually just before a Newton step).
     * It can be used to update some coefficients that will be used on material points (for example, compute
     * the parameters mu and lambda from the young modulus and poisson ratio given as data arguments).
     */
    void before_update() override {
        c01 = d_c01.getValue();
        c10 = d_c10.getValue();
        k   = d_k.getValue();
    }

    /**
     * Get the strain energy density Psi from the right Cauchy-Green strain tensor C.
     *
     * I1 = trace(C)
     * I2 = (trace(C))^2 - trace(C^2)
     * Psi(C) = C01 (I1 -3) + C10 * (I2 - 3) + K/2 (ln(J))^2
     *
     */
    Real
    strain_energy_density(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension>  & C) const override {
        const auto I1 = C.trace();
        const auto C_square = C*C;
        const auto I2 = (C.trace()*C.trace() - C_square.trace())/2;
        const auto lnJ = log(J);
        return c01 * (I1 - 3) + c10 * (I2 - 3) + k/2 * lnJ * lnJ;
        
    }

    /** Get the second Piola-Kirchhoff stress tensor from the right Cauchy-Green strain tensor C. */
    virtual Eigen::Matrix<Real, Dimension, Dimension>
    PK2_stress(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension>  & C) const override {
        const auto I1 = C.trace();
        static const auto Id = Eigen::Matrix<Real, Dimension, Dimension, Eigen::RowMajor>::Identity();
        const auto Ci = C.inverse();
        return 2*(c01 + c10*I1)*Id - 2*c10*C + Ci*(k*log(J));
    }

    /** Get the jacobian of the second Piola-Kirchhoff stress tensor w.r.t the right Cauchy-Green strain tensor C. */
    virtual Eigen::Matrix<Real, 6, 6>
    PK2_stress_jacobian(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension> & C) const override {
        using caribou::algebra::symmetric_dyad_1;
        using caribou::algebra::symmetric_dyad_2;

        static const auto Id = Eigen::Matrix<Real, Dimension, Dimension, Eigen::RowMajor>::Identity();
        const auto C_square = C*C;
        const auto C_square_inverse = C_square.inverse();
        const auto Ci = C.inverse().eval();
        return k*Ci*((1/2) * Id - Ci*log(J));

    }

private:
    // Private members
    Real c01; // Lame's mu parameter
    Real c10;  // Lame's lambda parameter
    Real k;

    // Data members
    sofa::core::objectmodel::Data<Real> d_c01;
    sofa::core::objectmodel::Data<Real> d_c10;
    sofa::core::objectmodel::Data<Real> d_k;
};

} // namespace SofaCaribou::material
