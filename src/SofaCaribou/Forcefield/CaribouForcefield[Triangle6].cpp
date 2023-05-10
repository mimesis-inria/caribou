#include <SofaCaribou/Forcefield/CaribouForcefield[Triangle6].h>
#include <SofaCaribou/Forcefield/CaribouForcefield.inl>
#include <SofaCaribou/Topology/CaribouTopology[Triangle6].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

using sofa::core::RegisterObject;
using namespace sofa::core::objectmodel;
using namespace caribou::geometry;
using namespace caribou;

namespace SofaCaribou::forcefield {

// -----------------------------------
// Triangle quadratic specialization
// -----------------------------------

// 2D

template <>
auto CaribouForcefield<Triangle6<_2D>>::templateName(const CaribouForcefield<Triangle6<_2D>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Triangle6<_2D>>::templateName();
}

template <>
void CaribouForcefield<Triangle6<_2D>>::triangulate_face(const Triangle6<_2D> & e, const std::size_t & /*face_id*/, std::vector<sofa::type::Vector3> & triangles_nodes) {
    auto triangle = std::vector {e.node(0), e.node(1), e.node(2)};
    for (const auto &n : triangle) {
        triangles_nodes.emplace_back(n[0], n[1], 0);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Triangle6<_2D>>;

// 3D

template <>
auto CaribouForcefield<Triangle6<_3D>>::templateName(const CaribouForcefield<Triangle6<_3D>> *) -> std::string {
    return SofaCaribou::topology::CaribouTopology<Triangle6<_3D>>::templateName();
}

template <>
void CaribouForcefield<Triangle6<_3D>>::triangulate_face(const Triangle6<_3D> & e, const std::size_t & /*face_id*/, std::vector<sofa::type::Vec3> & triangles_nodes) {
    auto triangle = std::vector {e.node(0), e.node(1), e.node(2)};
    for (const auto &n : triangle) {
        triangles_nodes.emplace_back(n[0], n[1], n[2]);
    }
}

// This will force the compiler to compile the following templated class
template class CaribouForcefield<Triangle6<_3D>>;

} // namespace SofaCaribou::forcefield
