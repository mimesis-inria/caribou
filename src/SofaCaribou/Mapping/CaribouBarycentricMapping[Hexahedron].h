#pragma once

#include <Caribou/Geometry/Hexahedron.h>

#include <SofaCaribou/Mapping/CaribouBarycentricMapping.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/defaulttype/VecTypes.h>
DISABLE_ALL_WARNINGS_END



namespace SofaCaribou::mapping {

/* using DataTypes = typename SofaCaribou::topology::SofaVecType<caribou::geometry::traits<caribou::geometry::Hexahedron<caribou::Linear>>::Dimension>::Type;
using MappedDataVecCoord = sofa::core::objectmodel::Data<sofa::defaulttype::Vec3Types::VecCoord>;
using DataVecCoord = sofa::core::objectmodel::Data<typename DataTypes::VecCoord>;
 */
// Hexahedron linear specialization
/* template <> void CaribouBarycentricMapping<caribou::geometry::Hexahedron<caribou::Linear>, sofa::defaulttype::Vec3Types>::apply(
    const sofa::core::MechanicalParams * /*mparams*,
    MappedDataVecCoord & data_output_mapped_position,
    const DataVecCoord & data_input_position); */
extern template class CaribouBarycentricMapping<caribou::geometry::Hexahedron<caribou::Linear>, sofa::defaulttype::Vec3Types>;
extern template class CaribouBarycentricMapping<caribou::geometry::Hexahedron<caribou::Linear>, sofa::defaulttype::Rigid3Types>;

// Hexahedron quadratic specialization
extern template class CaribouBarycentricMapping<caribou::geometry::Hexahedron<caribou::Quadratic>, sofa::defaulttype::Vec3Types>;
extern template class CaribouBarycentricMapping<caribou::geometry::Hexahedron<caribou::Quadratic>, sofa::defaulttype::Rigid3Types>;

} // namespace SofaCaribou::mapping