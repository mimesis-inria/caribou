#pragma once

#include <Caribou/constants.h>
#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/Mapping.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::mapping {

/**
 * Generic barycentric mapping.
 *
 * The CaribouBarycentricMapping allows to embed nodes into a container domain. For each embedded nodes, called
 * mapped nodes, the index of the element (from the container domain) that contains it is stored, paired to the
 * barycentric coordinates of the node within this element.
 *
 * @tparam Element The Element type of the container domain
 * @tparam MappedDataTypes The datatype of the embedded (mapped) degrees of freedom.
 */
template<typename Element, typename MappedDataTypes>
class CaribouBarycentricMapping
        : public sofa::core::Mapping<typename SofaCaribou::topology::SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type, MappedDataTypes> {
public:
    using DataTypes = typename SofaCaribou::topology::SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type;
    SOFA_CLASS(SOFA_TEMPLATE2(CaribouBarycentricMapping, Element, MappedDataTypes),
               SOFA_TEMPLATE2(sofa::core::Mapping, DataTypes, MappedDataTypes));

    // Type aliases
    using Domain = typename SofaCaribou::topology::CaribouTopology<Element>::Domain;

    static constexpr auto MappedDimension = MappedDataTypes::Coord::spatial_dimensions;
    using MappedScalar = typename MappedDataTypes::Real;
    using MappedDataVecCoord = sofa::core::objectmodel::Data<typename MappedDataTypes::VecCoord>;
    using MappedDataVecDeriv = sofa::core::objectmodel::Data<typename MappedDataTypes::VecDeriv>;
    using MappedDataMapMapSparseMatrix = sofa::core::objectmodel::Data<typename MappedDataTypes::MatrixDeriv>;

    static constexpr auto Dimension = DataTypes::Coord::spatial_dimensions;
    using Scalar = typename DataTypes::Real;
    using DataVecCoord = sofa::core::objectmodel::Data<typename DataTypes::VecCoord>;
    using DataVecDeriv = sofa::core::objectmodel::Data<typename DataTypes::VecDeriv>;
    using DataMapMapSparseMatrix  = sofa::core::objectmodel::Data<typename DataTypes::MatrixDeriv>;

    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<CaribouBarycentricMapping<Element, MappedDataTypes>, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    // Constructor
    CaribouBarycentricMapping();

    // Public virtual methods
     void init    () override;
     void apply   (const sofa::core::MechanicalParams* mparams, MappedDataVecCoord & output_mapped_position, const DataVecCoord& input_position) override;
     void applyJ  (const sofa::core::MechanicalParams* mparams, MappedDataVecDeriv & output_mapped_velocity, const DataVecDeriv& input_velocity) override;
     void applyJT (const sofa::core::MechanicalParams* mparams, DataVecDeriv & output_force, const MappedDataVecDeriv & input_mapped_force) override;
     void applyJT (const sofa::core::ConstraintParams* cparams, DataMapMapSparseMatrix & output_jacobian, const MappedDataMapMapSparseMatrix & input_jacobian) override;
     void draw    (const sofa::core::visual::VisualParams* vparams) override;

    static auto
    GetCustomTemplateName() -> std::string {
        return templateName();
    }

    static auto templateName(const CaribouBarycentricMapping<Element, MappedDataTypes>* = nullptr) -> std::string {
        return SofaCaribou::topology::CaribouTopology<Element>::templateName();
    }

    template <typename Derived>
    
    static auto canCreate(Derived * o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) -> bool;

private:
    // Data members
    Link<SofaCaribou::topology::CaribouTopology<Element>> d_topology;

    // Private members
    std::unique_ptr<caribou::topology::BarycentricContainer<Domain>> p_barycentric_container;

    Eigen::SparseMatrix<Scalar> p_J; ///< Mapping matrix
};

} // namespace SofaCaribou::mapping
