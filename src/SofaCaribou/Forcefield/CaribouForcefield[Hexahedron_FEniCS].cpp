#include <SofaCaribou/Forcefield/CaribouForcefield[Hexahedron_FEniCS].h>
#include <SofaCaribou/Forcefield/CaribouForcefield.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::objectmodel;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// --------------------------------
// Hexahedron_FEniCS linear specialization
// --------------------------------

template <>
auto CaribouForcefield<Hexahedron_FEniCS<Linear>>::templateName(const CaribouForcefield<Hexahedron_FEniCS<Linear>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Hexahedron_FEniCS<Linear>>::templateName();
}

template <>
void CaribouForcefield<Hexahedron_FEniCS<Linear>>::triangulate_face(const Hexahedron_FEniCS < Linear> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes) {
    const auto face_indices = e.boundary_elements_node_indices();
    auto triangle_1 = std::vector {e.node(face_indices[face_id][0]), e.node(face_indices[face_id][1]), e.node(face_indices[face_id][2])};
    auto triangle_2 = std::vector {e.node(face_indices[face_id][2]), e.node(face_indices[face_id][3]), e.node(face_indices[face_id][0])};
    for (const auto &n : triangle_1) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
    for (const auto &n : triangle_2) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Hexahedron_FEniCS < Linear>>;

// -----------------------------------
// Hexahedron_FEniCS quadratic specialization
// -----------------------------------

template <>
auto CaribouForcefield<Hexahedron_FEniCS<Quadratic>>::templateName(const CaribouForcefield<Hexahedron_FEniCS<Quadratic>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Hexahedron_FEniCS<Quadratic>>::templateName();
}

template <>
void CaribouForcefield<Hexahedron_FEniCS<Quadratic>>::triangulate_face(const Hexahedron_FEniCS < Quadratic> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes) {
    const auto face_indices = e.boundary_elements_node_indices();
    auto triangle_1 = std::vector {e.node(face_indices[face_id][0]), e.node(face_indices[face_id][1]), e.node(face_indices[face_id][2])};
    auto triangle_2 = std::vector {e.node(face_indices[face_id][2]), e.node(face_indices[face_id][3]), e.node(face_indices[face_id][0])};
    for (const auto &n : triangle_1) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
    for (const auto &n : triangle_2) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Hexahedron_FEniCS < Quadratic>>;

} // namespace SofaCaribou::forcefield
