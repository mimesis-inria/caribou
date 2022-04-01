#include <SofaCaribou/Forcefield/CaribouForcefield[Quad8].h>
#include <SofaCaribou/Forcefield/CaribouForcefield.inl>
#include <SofaCaribou/Topology/CaribouTopology[Quad8].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::objectmodel;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// -----------------------------------
// Quad quadratic specialization
// -----------------------------------

// 2D

template <>
auto CaribouForcefield<Quad8<_2D>>::templateName(const CaribouForcefield<Quad8<_2D>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Quad8<_2D>>::templateName();
}

template <>
void CaribouForcefield<Quad8<_2D>>::triangulate_face(const Quad8<_2D> & e, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & triangles_nodes) {
    auto triangle_1 = std::vector {e.node(0), e.node(1), e.node(2)};
    auto triangle_2 = std::vector {e.node(2), e.node(3), e.node(0)};
    for (const auto &n : triangle_1) {
        triangles_nodes.emplace_back(n[0], n[1], 0);
    }
    for (const auto &n : triangle_2) {
        triangles_nodes.emplace_back(n[0], n[1], 0);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Quad8<_2D>>;

// 3D

template <>
auto CaribouForcefield<Quad8<_3D>>::templateName(const CaribouForcefield<Quad8<_3D>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Quad8<_3D>>::templateName();
}

template <>
void CaribouForcefield<Quad8<_3D>>::triangulate_face(const Quad8<_3D> & e, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & triangles_nodes) {
    auto triangle_1 = std::vector {e.node(0), e.node(1), e.node(2)};
    auto triangle_2 = std::vector {e.node(2), e.node(3), e.node(0)};
    for (const auto &n : triangle_1) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
    for (const auto &n : triangle_2) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Quad8<_3D>>;

} // namespace SofaCaribou::forcefield
