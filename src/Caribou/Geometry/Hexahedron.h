#ifndef CARIBOU_GEOMETRY_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_HEXAHEDRON_H

#include <Caribou/Geometry/Entity.h>

namespace caribou
{
namespace geometry
{

/* Linear hexahedron node numbering
 *    7-------6
 *   /|      /|
 *  / |     / |
 * 3--|----2  |
 * |  4----|--5
 * | /     | /
 * 0-------1
 */

struct LinearHexahedron : public Entity
{

    std::array<PointType, 8> nodes;
};

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_HEXAHEDRON_H
