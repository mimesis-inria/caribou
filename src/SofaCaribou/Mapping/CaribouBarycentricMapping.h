#pragma once

#include <Caribou/constants.h>
#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.h>
#include <Caribou/Mechanics/Elasticity/Strain.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/Mapping.h>
#include <sofa/helper/OptionsGroup.h>
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
    static constexpr auto MapTotalSize = MappedDataTypes::Coord::total_size;
    static constexpr auto MapForce = 6;
    
    using MappedScalar = typename MappedDataTypes::Real;
    using MappedDataVecCoord = sofa::core::objectmodel::Data<typename MappedDataTypes::VecCoord>;
    using MappedDataVecDeriv = sofa::core::objectmodel::Data<typename MappedDataTypes::VecDeriv>;
    using MappedDataMapMapSparseMatrix = sofa::core::objectmodel::Data<typename MappedDataTypes::MatrixDeriv>;

    using Mat3x3   = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
    using Vect3x1   = Eigen::Vector3d;
    static constexpr INTEGER_TYPE NumberOfNodesPerElement = caribou::geometry::traits<Element>::NumberOfNodesAtCompileTime;

    using LocalCoordinates = typename Element::LocalCoordinates;

    static constexpr auto Dimension = DataTypes::Coord::spatial_dimensions;
    using Scalar = typename DataTypes::Real;
    using DataVecCoord = sofa::core::objectmodel::Data<typename DataTypes::VecCoord>;
    using DataVecDeriv = sofa::core::objectmodel::Data<typename DataTypes::VecDeriv>;
    using DataMapMapSparseMatrix  = sofa::core::objectmodel::Data<typename DataTypes::MatrixDeriv>;

    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<CaribouBarycentricMapping<Element, MappedDataTypes>, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    template <typename T>
    using Data = sofa::core::objectmodel::Data<T>;

    /// Rotation extraction methods
    enum RotationExtractionMethod {
        /// An orthogonal frame is built by projecting the mapping point onto the element's faces.
        Frame = 0,

        /// The rotation is extracted from the deformation tensor evaluated at the mapped point
        /// using the Singular Value Decomposition (SVD) method.
        SVD = 1,

        /// The rotation is extracted from the deformation tensor using the Analytic Polar
        /// Decomposition (APD) method described in "Fast Corotated FEM using Operator Splitting"
        /// T. Kugelstadt et al. 2018.
        APD = 2,
    };

    // Constructor
    CaribouBarycentricMapping();

    // Public virtual methods
    CARIBOU_API void init    () override;
    CARIBOU_API void apply   (const sofa::core::MechanicalParams* mparams, MappedDataVecCoord & output_mapped_position, const DataVecCoord& input_position) override;
    CARIBOU_API void applyJ  (const sofa::core::MechanicalParams* mparams, MappedDataVecDeriv & output_mapped_velocity, const DataVecDeriv& input_velocity) override;
    CARIBOU_API void applyJT (const sofa::core::MechanicalParams* mparams, DataVecDeriv & output_force, const MappedDataVecDeriv & input_mapped_force) override;
    CARIBOU_API void applyJT (const sofa::core::ConstraintParams* cparams, DataMapMapSparseMatrix & output_jacobian, const MappedDataMapMapSparseMatrix & input_jacobian) override;
    CARIBOU_API void draw    (const sofa::core::visual::VisualParams* vparams) override;

    [[nodiscard]] CARIBOU_API auto
    getTemplateName() const -> std::string override {
        return templateName(this);
    }

    static auto templateName(const CaribouBarycentricMapping<Element, MappedDataTypes>* = nullptr) -> std::string {
        return SofaCaribou::topology::CaribouTopology<Element>::templateName() + "," + MappedDataTypes::Name();
    }

    template <typename Derived>
    CARIBOU_API
    static auto canCreate(Derived * o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) -> bool;

private:
    // Data members
    Link<SofaCaribou::topology::CaribouTopology<Element>> d_topology;

    Data< sofa::helper::OptionsGroup > d_rotation_extraction_method;

    // Private members
    std::unique_ptr<caribou::topology::BarycentricContainer<Domain>> p_barycentric_container;

    Eigen::SparseMatrix<Scalar> p_J; ///< Mapping matrix

    std::vector<Mat3x3> initial_intersections;

};

} // namespace SofaCaribou::mapping