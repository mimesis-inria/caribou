#include <SofaCaribou/Forcefield/HyperelasticForcefield[Hexahedron].h>
#include <SofaCaribou/Forcefield/HyperelasticForcefield.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::topology;
using namespace sofa::component::topology;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// --------------------------------
// Hexahedron linear specialization
// --------------------------------

template <>
void HyperelasticForcefield<Hexahedron<Linear>>::create_domain_from(BaseMeshTopology * topology) {
    // A mesh is required in order to create the domain
    if (not p_mesh) {
        return;
    }

    // Get a reference to the list of hexahedrons
    const auto & hexahedrons = topology->getHexahedra();
    const auto * data = hexahedrons.data()->data();

    // Create the Domain instance
    using PointID = sofa::core::topology::Topology::PointID;
    p_domain = p_mesh->add_domain<Hexahedron<Linear>, PointID> (
        data, hexahedrons.size(), Hexahedron<Linear>::NumberOfNodesAtCompileTime
    );
}

template <>
auto HyperelasticForcefield<Hexahedron<Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
        dynamic_cast<const HexahedronSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto HyperelasticForcefield<Hexahedron<Linear>>::templateName(const HyperelasticForcefield<Hexahedron<Linear>> *) -> std::string {
    return "Hexahedron";
}

// This will force the compiler to compile the following templated class
template class HyperelasticForcefield<Hexahedron < Linear>>;

// -----------------------------------
// Hexahedron quadratic specialization
// -----------------------------------
template <>
void HyperelasticForcefield<Hexahedron<Quadratic>>::create_domain_from(BaseMeshTopology *) {
    // We cannot create a quadratic topology from a linear one
}

template <>
auto HyperelasticForcefield<Hexahedron<Quadratic>>::mesh_is_compatible(const BaseMeshTopology *) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

template <>
auto HyperelasticForcefield<Hexahedron<Quadratic>>::templateName(const HyperelasticForcefield<Hexahedron<Quadratic>> *) -> std::string {
    return "Hexahedron20";
}

template <>
void HyperelasticForcefield<Hexahedron<Quadratic>>::draw(const sofa::core::visual::VisualParams * vparams) {
    using Color = sofa::helper::types::RGBAColor;
    constexpr static auto NumberOfFaces = 6;

    if (!vparams->displayFlags().getShowForceFields())
        return;

    const auto nb_elements = number_of_elements();

    if (nb_elements == 0)
        return;

    vparams->drawTool()->saveLastState();
    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,true);
    vparams->drawTool()->disableLighting();

    const double & scale = d_drawScale.getValue();
    const VecCoord& sofa_x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.data()->data(),  sofa_x.size(), Dimension);

    std::vector< sofa::defaulttype::Vec<Dimension, Real> > faces_points[NumberOfFaces];
    for (std::size_t face_id = 0; face_id < NumberOfFaces; ++face_id) {
        faces_points[face_id].reserve(nb_elements*4);
    }

    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        // Fetch the node indices of the element
        auto node_indices = p_domain->element_indices(element_id);
        Matrix<NumberOfNodes, Dimension> element_nodes_position;

        // Fetch the initial positions of the element's nodes
        for (std::size_t node_id = 0; node_id < NumberOfNodes; ++node_id) {
            element_nodes_position.row(node_id) = X.row(node_indices[node_id]);
        }

        // Create an Element instance from the node positions
        const auto e = Hexahedron<Quadratic>(element_nodes_position);

        // Scale down the element around its center point
        const auto c = e.center();
        for (std::size_t node_id = 0; node_id < NumberOfNodes; ++node_id) {
            const auto & p = element_nodes_position.row(node_id).transpose();
            element_nodes_position.row(node_id) = (c + (p - c)*scale).transpose();
        }

        // Push the faces scaled-down nodes
        const auto face_node_indices = e.boundary_elements_node_indices();
        for (std::size_t face_id = 0; face_id < NumberOfFaces; ++face_id) {
            for (std::size_t face_node_id = 0; face_node_id < 4; ++face_node_id) {
                const auto & p = element_nodes_position.row(face_node_indices[face_id][face_node_id]);
                faces_points[face_id].emplace_back(p[0], p[1], p[2]);
            }
        }

        for (std::size_t face_id = 0; face_id < NumberOfFaces; ++face_id) {
            const auto &hex_color = kelly_colors_hex[face_id % 20];
            const Color face_color(
                    static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(16))) / 255.),
                    static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(8))) / 255.),
                    static_cast<float> ((static_cast<unsigned char> (hex_color >> static_cast<unsigned>(0))) / 255.),
                    static_cast<float> (1));

            vparams->drawTool()->drawQuads(faces_points[face_id], face_color);
        }
    }

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0,false);

    vparams->drawTool()->restoreLastState();
}

// This will force the compiler to compile the following templated class
template class HyperelasticForcefield<Hexahedron < Quadratic>>;

} // namespace SofaCaribou::forcefield

namespace sofa::core::objectmodel {
using namespace SofaCaribou::forcefield;

[[maybe_unused]]
static int _c_ = RegisterObject("Caribou hyperelastic force field")
    .add<HyperelasticForcefield<Hexahedron<Linear>>>()
    .add<HyperelasticForcefield<Hexahedron<Quadratic>>>();
}
