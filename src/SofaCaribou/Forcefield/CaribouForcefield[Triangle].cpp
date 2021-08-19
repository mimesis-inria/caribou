#include <SofaCaribou/Forcefield/CaribouForcefield[Triangle].h>
#include <SofaCaribou/Forcefield/CaribouForcefield.inl>
#include <SofaCaribou/Topology/CaribouTopology[Triangle].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::objectmodel;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// --------------------------------
// Triangle linear specialization
// --------------------------------

// 2D

template <>
auto CaribouForcefield<Triangle<_2D, Linear>>::templateName(const CaribouForcefield<Triangle<_2D, Linear>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Triangle<_2D, Linear>>::templateName();
}

template <>
void CaribouForcefield<Triangle<_2D, Linear>>::triangulate_face(const Triangle<_2D, Linear> & e, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & triangles_nodes) {
    auto triangle = std::vector {e.node(0), e.node(1), e.node(2)};
    for (const auto &n : triangle) {
        triangles_nodes.emplace_back(n[0], n[1], 0);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Triangle<_2D, Linear>>;

// 3D

template <>
auto CaribouForcefield<Triangle<_3D, Linear>>::templateName(const CaribouForcefield<Triangle<_3D, Linear>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Triangle<_3D, Linear>>::templateName();
}

template <>
void CaribouForcefield<Triangle<_3D, Linear>>::triangulate_face(const Triangle<_3D, Linear> & e, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & triangles_nodes) {
    auto triangle = std::vector {e.node(0), e.node(1), e.node(2)};
    for (const auto &n : triangle) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Triangle<_3D, Linear>>;

// -----------------------------------
// Triangle quadratic specialization
// -----------------------------------

// 2D

template <>
auto CaribouForcefield<Triangle<_2D, Quadratic>>::templateName(const CaribouForcefield<Triangle<_2D, Quadratic>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Triangle<_2D, Quadratic>>::templateName();
}

template <>
void CaribouForcefield<Triangle<_2D, Quadratic>>::triangulate_face(const Triangle<_2D, Quadratic> & e, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & triangles_nodes) {
    auto triangle = std::vector {e.node(0), e.node(1), e.node(2)};
    for (const auto &n : triangle) {
        triangles_nodes.emplace_back(n[0], n[1], 0);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Triangle<_2D, Quadratic>>;

// 3D

template <>
auto CaribouForcefield<Triangle<_3D, Quadratic>>::templateName(const CaribouForcefield<Triangle<_3D, Quadratic>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Triangle<_3D, Quadratic>>::templateName();
}

template <>
void CaribouForcefield<Triangle<_3D, Quadratic>>::triangulate_face(const Triangle<_3D, Quadratic> & e, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & triangles_nodes) {
    auto triangle = std::vector {e.node(0), e.node(1), e.node(2)};
    for (const auto &n : triangle) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Triangle<_3D, Quadratic>>;

} // namespace SofaCaribou::forcefield
