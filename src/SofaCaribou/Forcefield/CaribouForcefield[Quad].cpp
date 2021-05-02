#include <SofaCaribou/Forcefield/CaribouForcefield[Quad].h>
#include <SofaCaribou/Forcefield/CaribouForcefield.inl>
#include <SofaCaribou/Topology/CaribouTopology[Quad].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseTopology/QuadSetTopologyContainer.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::topology;
using namespace sofa::core::objectmodel;
using namespace sofa::component::topology;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// --------------------------------
// Quad linear specialization
// --------------------------------

// 2D

template <>
auto CaribouForcefield<Quad<_2D, Linear>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> BaseData * {
    return topology->findData("quads");
}

template <>
auto CaribouForcefield<Quad<_2D, Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const QuadSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouForcefield<Quad<_2D, Linear>>::templateName(const CaribouForcefield<Quad<_2D, Linear>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Quad<_2D, Linear>>::templateName();
}

template <>
void CaribouForcefield<Quad<_2D, Linear>>::triangulate_face(const Quad<_2D, Linear> & e, const std::size_t & /*face_id*/, std::vector<sofa::defaulttype::Vector3> & triangles_nodes) {
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
template class CaribouForcefield<Quad<_2D, Linear>>;

// 3D

template <>
auto CaribouForcefield<Quad<_3D, Linear>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * topology) -> BaseData * {
    return topology->findData("quads");
}

template <>
auto CaribouForcefield<Quad<_3D, Linear>>::mesh_is_compatible(const BaseMeshTopology * topology) -> bool {
    return (
            dynamic_cast<const QuadSetTopologyContainer*>(topology) != nullptr
    );
}

template <>
auto CaribouForcefield<Quad<_3D, Linear>>::templateName(const CaribouForcefield<Quad<_3D, Linear>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Quad<_3D, Linear>>::templateName();
}

template <>
void CaribouForcefield<Quad<_3D, Linear>>::triangulate_face(const Quad<_3D, Linear> & e, const std::size_t & /*face_id*/, std::vector<sofa::defaulttype::Vector3> & triangles_nodes) {
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
template class CaribouForcefield<Quad<_3D, Linear>>;

// -----------------------------------
// Quad quadratic specialization
// -----------------------------------

// 2D

template <>
auto CaribouForcefield<Quad<_2D, Quadratic>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * /*topology*/) -> BaseData * {
    // We cannot create a quadratic topology from a linear one
    return nullptr;
}

template <>
auto CaribouForcefield<Quad<_2D, Quadratic>>::mesh_is_compatible(const BaseMeshTopology * /*topology*/) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

template <>
auto CaribouForcefield<Quad<_2D, Quadratic>>::templateName(const CaribouForcefield<Quad<_2D, Quadratic>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Quad<_2D, Quadratic>>::templateName();
}

template <>
void CaribouForcefield<Quad<_2D, Quadratic>>::triangulate_face(const Quad<_2D, Quadratic> & e, const std::size_t & /*face_id*/, std::vector<sofa::defaulttype::Vector3> & triangles_nodes) {
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
template class CaribouForcefield<Quad<_2D, Quadratic>>;

// 3D

template <>
auto CaribouForcefield<Quad<_3D, Quadratic>>::get_indices_from(const sofa::core::topology::BaseMeshTopology * /*topology*/) -> BaseData * {
    // We cannot create a quadratic topology from a linear one
    return nullptr;
}

template <>
auto CaribouForcefield<Quad<_3D, Quadratic>>::mesh_is_compatible(const BaseMeshTopology * /*topology*/) -> bool {
    // We cannot create a quadratic topology from a linear one
    return false;
}

template <>
auto CaribouForcefield<Quad<_3D, Quadratic>>::templateName(const CaribouForcefield<Quad<_3D, Quadratic>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Quad<_3D, Quadratic>>::templateName();
}

template <>
void CaribouForcefield<Quad<_3D, Quadratic>>::triangulate_face(const Quad<_3D, Quadratic> & e, const std::size_t & /*face_id*/, std::vector<sofa::defaulttype::Vector3> & triangles_nodes) {
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
template class CaribouForcefield<Quad<_3D, Quadratic>>;

} // namespace SofaCaribou::forcefield
