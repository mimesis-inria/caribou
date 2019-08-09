#ifndef CARIBOU_GEOMETRY_INTERNAL_BASETETRAHEDRON_H
#define CARIBOU_GEOMETRY_INTERNAL_BASETETRAHEDRON_H

#include <Caribou/config.h>

namespace caribou::geometry::internal {

template<typename CanonicalElementType, typename TetrahedronType>
struct BaseTetrahedron : public CanonicalElementType
{
    static constexpr INTEGER_TYPE NumberOfNodes = TetrahedronType::NumberOfNodes;
    using CanonicalElement = CanonicalElementType;
    using LocalCoordinates = typename CanonicalElement::LocalCoordinates;
    using WorldCoordinates = Eigen::Matrix<FLOATING_POINT_TYPE, 3, 1>;
};

} // namespace caribou::geometry::internal

#endif //CARIBOU_GEOMETRY_INTERNAL_BASETETRAHEDRON_H
