#pragma once

#include <Caribou/Algebra/Tensor.h>
#include <SofaCaribou/Material/HyperelasticMaterial.h>

namespace SofaCaribou::material {

template<class DataTypes>
class NeoHookeanMaterial : public HyperelasticMaterial<DataTypes> {
    static constexpr auto Dimension = DataTypes::spatial_dimensions;
    using Coord = typename DataTypes::Coord;
    using Real  = typename Coord::value_type;
public:
    SOFA_CLASS(SOFA_TEMPLATE(NeoHookeanMaterial, DataTypes), SOFA_TEMPLATE(HyperelasticMaterial, DataTypes));

    NeoHookeanMaterial()
        : d_young_modulus(initData(&d_young_modulus,
            Real(1000), "young_modulus",
            "Young's modulus of the material",
            true /*displayed_in_GUI*/, false /*read_only_in_GUI*/))
        , d_poisson_ratio(initData(&d_poisson_ratio,
            Real(0.3),  "poisson_ratio",
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
        const Real young_modulus = d_young_modulus.getValue();
        const Real poisson_ratio = d_poisson_ratio.getValue();
        mu = young_modulus / (2.0 * (1.0 + poisson_ratio));
        l = young_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    }

    /**
     * Get the strain energy density Psi from the right Cauchy-Green strain tensor C.
     *
     * Psi(C) = mu/2 (tr(C)-3) - mu ln(J) + lambda/2 (ln(J))^2
     *
     */
    Real
    strain_energy_density(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension>  & C) const override {
        const auto lnJ = log(J);
        return mu/2.*(C.trace()-3) - mu*lnJ + l/2 *lnJ*lnJ;
    }

    /** Get the second Piola-Kirchhoff stress tensor from the right Cauchy-Green strain tensor C. */
    virtual Eigen::Matrix<Real, Dimension, Dimension>
    PK2_stress(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension>  & C) const override {

        static const auto Id = Eigen::Matrix<Real, Dimension, Dimension, Eigen::RowMajor>::Identity();
        const auto Ci = C.inverse();

        return (Id - Ci)*mu + Ci*(l*log(J));
    }

    /** Get the jacobian of the second Piola-Kirchhoff stress tensor w.r.t the right Cauchy-Green strain tensor C. */
    virtual Eigen::Matrix<Real, 6, 6>
    PK2_stress_jacobian(const Real & J, const Eigen::Matrix<Real, Dimension, Dimension> & C) const override {
        using caribou::algebra::symmetric_dyad_1;
        using caribou::algebra::symmetric_dyad_2;

        const auto Ci = C.inverse().eval();

        Eigen::Matrix<Real, 6, 6> D = l*symmetric_dyad_1(Ci) + 2*(mu - l*log(J))*symmetric_dyad_2(Ci);
        return D;
    }

<<<<<<< HEAD

=======
>>>>>>> FenicCS-features

private:
    // Private members
    Real mu; // Lame's mu parameter
    Real l;  // Lame's lambda parameter

    // Data members
    sofa::core::objectmodel::Data<Real> d_young_modulus;
    sofa::core::objectmodel::Data<Real> d_poisson_ratio;
};

} // namespace SofaCaribou::material
