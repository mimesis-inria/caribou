#include <SofaCaribou/Forcefield/CaribouForcefield[Tetrahedron].h>
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
// Tetrahedron linear specialization
// --------------------------------
template <>
auto CaribouForcefield<Tetrahedron<Linear>>::templateName(const CaribouForcefield<Tetrahedron<Linear>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Tetrahedron<Linear>>::templateName();
}

template<>
void CaribouForcefield<Tetrahedron<Linear>>::triangulate_face(const Tetrahedron<Linear> &e, const std::size_t &face_id, std::vector<sofa::type::Vector3> &triangles_nodes) {
    const auto face_indices = e.boundary_elements_node_indices();
    auto triangle_1 = std::array{e.node(face_indices[face_id][0]), e.node(face_indices[face_id][1]), e.node(face_indices[face_id][2])};
    for (const auto &n : triangle_1) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Tetrahedron < Linear>>;

// -----------------------------------
// Tetrahedron quadratic specialization
// -----------------------------------

template <>
auto CaribouForcefield<Tetrahedron<Quadratic>>::templateName(const CaribouForcefield<Tetrahedron<Quadratic>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Tetrahedron<Quadratic>>::templateName();
}

template<>
void CaribouForcefield<Tetrahedron<Quadratic>>::triangulate_face(const Tetrahedron<Quadratic> &e, const std::size_t &face_id, std::vector<sofa::type::Vector3> &triangles_nodes) {
    const auto face_indices = e.boundary_elements_node_indices();
    auto triangle_1 = std::array{e.node(face_indices[face_id][0]), e.node(face_indices[face_id][1]), e.node(face_indices[face_id][2])};
    for (const auto &n : triangle_1) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Tetrahedron < Quadratic>>;

} // namespace SofaCaribou::forcefield
