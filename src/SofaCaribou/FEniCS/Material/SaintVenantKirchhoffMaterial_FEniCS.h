#pragma once

#include <SofaCaribou/FEniCS/Material/HyperelasticMaterial_FEniCS.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <SofaCaribou/FEniCS/Material/SaintVenantKirchhoff_Tetra.h>
#include <SofaCaribou/FEniCS/Material/SaintVenantKirchhoff_Tetra_Order2.h>
#include <SofaCaribou/FEniCS/Material/SaintVenantKirchhoff_Hexa.h>
#include <SofaCaribou/FEniCS/Material/SaintVenantKirchhoff_Hexa_Order2.h>
#include <iostream>


namespace SofaCaribou::material {

template<typename Element, typename DataTypes>
class SaintVenantKirchhoffMaterial_FEniCS : public HyperelasticMaterial_FEniCS<Element, DataTypes> {
    static constexpr auto Dimension = DataTypes::spatial_dimensions;
    static constexpr INTEGER_TYPE NumberOfNodesPerElement = caribou::geometry::traits<Element>::NumberOfNodesAtCompileTime;
    using Coord = typename DataTypes::Coord;
    using Real  = typename Coord::value_type;
public:
    SOFA_CLASS(SOFA_TEMPLATE2(SaintVenantKirchhoffMaterial_FEniCS, Element, DataTypes), SOFA_TEMPLATE2(HyperelasticMaterial_FEniCS, Element, DataTypes));

    SaintVenantKirchhoffMaterial_FEniCS()
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

        Eigen::Matrix<Real, 2*Dimension, 2*Dimension> D;
        D <<
              l + 2*mu,    l,          l,       0,  0,  0,
                l,       l + 2*mu,     l,       0,  0,  0,
                l,         l,        l + 2*mu,  0,  0,  0,
                0,         0,          0,      mu,  0,  0,
                0,         0,          0,       0, mu,  0,
                0,         0,          0,       0,  0, mu;
        C = D;
    }

    /**
     * Get the strain energy density Psi from the right Cauchy-Green strain tensor C.
     *
     * Psi(E) = lambda/2 (tr(E))^2 + mu tr(E*E)
     *
     */
    Real
    strain_energy_density(const Real & /*J*/, const Eigen::Matrix<Real, Dimension, Dimension>  & C) const override {
        static const auto Id = Eigen::Matrix<Real, Dimension, Dimension, Eigen::RowMajor>::Identity();
        const auto E = (1/2. * (C - Id)).eval();
        const auto trE  = E.trace();
        const auto trEE = (E*E).trace();
        return l/2.*(trE*trE) + mu*trEE;
    }

    /** Get the second Piola-Kirchhoff stress tensor from the right Cauchy-Green strain tensor C. */
    Eigen::Matrix<Real, Dimension, Dimension>
    PK2_stress(const Real & /*J*/, const Eigen::Matrix<Real, Dimension, Dimension>  & C) const override {
        static const auto Id = Eigen::Matrix<Real, Dimension, Dimension, Eigen::RowMajor>::Identity();
        const auto E = (1/2. * (C - Id)).eval();
        return l*E.trace()*Id + 2*mu*E;
    }

    /** Get the jacobian of the second Piola-Kirchhoff stress tensor w.r.t the right Cauchy-Green strain tensor C. */
    Eigen::Matrix<Real, 6, 6>
    PK2_stress_jacobian(const Real & /*J*/, const Eigen::Matrix<Real, Dimension, Dimension> & /*C*/) const override {
        return C;
    }


    void
    calculate_nodal_forces(Eigen::Matrix<Real, NumberOfNodesPerElement, Dimension> nodal_forces, Eigen::Matrix<Real, NumberOfNodesPerElement, Dimension> coefficients, Eigen::Matrix<Real, NumberOfNodesPerElement, Dimension> current_nodes_position) const override{
//    std::cout << "nodal forces" << nodal_forces << "\n";
//    std::cout << "coefficients" << coefficients << "\n";
//    std::cout << "position" << current_nodes_position << "\n";
//    integral->tabulate_tensor(nodal_forces.data(), coefficients.data(), constants, current_nodes_position.data(), nullptr, nullptr);
//    std::cout << "final nodal forces" << nodal_forces << "\n";
      std::cout << "Teeeeeeeeeeeeeeeeeeeeeeeeeeest";
      std::cout << NumberOfNodesPerElement << "/n";
    }



private:
    // Private members
    Real mu; // Lame's mu parameter
    Real l;  // Lame's lambda parameter
    Eigen::Matrix<Real, 2*Dimension, 2*Dimension> C; // Constant elasticity matrix (jacobian of the stress tensor)

    // Data members
    sofa::core::objectmodel::Data<Real> d_young_modulus;
    sofa::core::objectmodel::Data<Real> d_poisson_ratio;
//    const ufc_scalar_t constants[2] = {d_young_modulus.getValue(), d_poisson_ratio.getValue()};
//    // Get the single cell integral
//    const ufc_integral *integral =
//        form_SaintVenantKirchhoff_Tetra_F->integrals(ufc_integral_type::cell)[0];
};

} // namespace SofaCaribou::material
