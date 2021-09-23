#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END

#include <Caribou/config.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Caribou/Topology/Mesh.h>

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa { using Index = unsigned int; }
#endif

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210600)
namespace sofa::type { using Vector3 = ::sofa::defaulttype::Vector3; }
#endif

namespace SofaCaribou::forcefield {

// Traits to get the Sofa vector type from the dimension
template <std::size_t Dim> struct SofaVecType {};
template <> struct SofaVecType<1> { using Type = sofa::defaulttype::Vec1Types; };
template <> struct SofaVecType<2> { using Type = sofa::defaulttype::Vec2Types; };
template <> struct SofaVecType<3> { using Type = sofa::defaulttype::Vec3Types; };

template <typename Element>
class CaribouForcefield : public sofa::core::behavior::ForceField<typename SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(CaribouForcefield, Element), SOFA_TEMPLATE(sofa::core::behavior::ForceField, typename SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type));

    // Type aliases
    using DataTypes = typename SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using Real     = typename DataTypes::Real;

    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<CaribouForcefield<Element>, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    template<int nRows, int nColumns>
    using Matrix =Eigen::Matrix<Real, nRows, nColumns, Eigen::RowMajor>;

    template<int nRows>
    using Vector = Eigen::Matrix<Real, nRows, 1>;

    // Constant properties
    static constexpr INTEGER_TYPE Dimension = caribou::geometry::traits<Element>::Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesPerElement = caribou::geometry::traits<Element>::NumberOfNodesAtCompileTime;
    static constexpr INTEGER_TYPE NumberOfGaussNodesPerElement = caribou::geometry::traits<Element>::NumberOfGaussNodesAtCompileTime;

    // Public methods
    
    CaribouForcefield();

    
    void init() override;

    template <typename Derived>
    
    static auto canCreate(Derived * o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) -> bool;

    [[nodiscard]] auto
    getTemplateName() const -> std::string override {
        return templateName(this);
    }

    static auto templateName(const CaribouForcefield<Element>* = nullptr) -> std::string;

    
    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;

    
    void draw(const sofa::core::visual::VisualParams* vparams) override;

    /** Get the number of elements contained in this field **/
    [[nodiscard]] inline
    auto number_of_elements() const noexcept -> std::size_t {
        return (p_topology ? p_topology->number_of_elements() : 0);
    }

    [[nodiscard]] inline
    auto topology() const noexcept -> typename SofaCaribou::topology::CaribouTopology<Element>::SPtr {
        return p_topology;
    }

private:

    /**
     * Construct triangles representing the visual of a given face to be draw by the force field.
     * @param e [input] The element being drawn by the force field.
     * @param face_id [input] The id of the face being drawn by the force field.
     * @param triangles_nodes [output] The position of each triangle nodes that will have to be draw to correctly represent the
     *                        face (visually).
     */
    void triangulate_face(const Element & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);

    // Data members
    /// This link is specifically set to point towards a very general BaseObject since it can be either a
    /// BaseMeshTopology (SOFA topology container), or a CaribouTopology.
    Link<sofa::core::objectmodel::BaseObject> d_topology_container;
    sofa::core::objectmodel::Data<double> d_drawScale;

    // Private variables
    /// Pointer to a CaribouTopology. This pointer will be null if a CaribouTopology
    /// is found within the scene graph and linked using the d_topology_container data
    /// parameter. Otherwise, if a compatible SOFA's topology (see SofaCaribou::topology::CaribouTopology::mesh_is_compatible())
    /// is found and linked, an internal CaribouTopology component will be created
    /// and its pointer will be stored here.
    typename SofaCaribou::topology::CaribouTopology<Element>::SPtr p_topology;
};

} // namespace SofaCaribou::forcefield