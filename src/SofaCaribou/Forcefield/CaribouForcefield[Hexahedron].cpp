#include <SofaCaribou/Forcefield/CaribouForcefield[Hexahedron].h>
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
// Hexahedron linear specialization
// --------------------------------

template <>
auto CaribouForcefield<Hexahedron<Linear>>::templateName(const CaribouForcefield<Hexahedron<Linear>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Hexahedron<Linear>>::templateName();
}

template <>
void CaribouForcefield<Hexahedron<Linear>>::triangulate_face(const Hexahedron < Linear> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes) {
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
template class CaribouForcefield<Hexahedron < Linear>>;

// -----------------------------------
// Hexahedron quadratic specialization
// -----------------------------------

template <>
auto CaribouForcefield<Hexahedron<Quadratic>>::templateName(const CaribouForcefield<Hexahedron<Quadratic>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Hexahedron<Quadratic>>::templateName();
}

template <>
void CaribouForcefield<Hexahedron<Quadratic>>::triangulate_face(const Hexahedron < Quadratic> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes) {
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
template class CaribouForcefield<Hexahedron < Quadratic>>;

} // namespace SofaCaribou::forcefield
