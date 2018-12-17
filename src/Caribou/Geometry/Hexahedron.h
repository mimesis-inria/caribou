#ifndef CARIBOU_GEOMETRY_HEXAHEDRON_H
#define CARIBOU_GEOMETRY_HEXAHEDRON_H

#include <array>

namespace caribou
{
namespace geometry
{

/**
 * The base hexahedron class (of n nodes) in 3D space.
 *
 * This class will be later defined by an explicit hexahedron class (8-nodes linear, 20-nodes quadratic, etc.) and the
 * type (regular or non-regular).
 *
 * ** Do not use this class directly. Use instead caribou::geometry::*Hexahdron. **
 *
 * The functions declared in this class can be used with any type of hexahedron (linear, quadratic, regular and non-regular).
 *
 * Do to so, it uses the Curiously Recurring Template Pattern (CRTP) :
 *    https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
 */
template <size_t NNodes, typename HexahedronExplicitType_>
struct BaseHexahedron
{

    static constexpr size_t NumberOfNodes = NNodes;

    static_assert(NumberOfNodes >= 8, "A hexahedron must have at least eight nodes.");

    /** Delete the default constructor as it should be defined by an explicit hexahedron class **/
    Hexahedron () = delete;
};

} // namespace geometry

} // namespace caribou

#endif //CARIBOU_GEOMETRY_HEXAHEDRON_H
