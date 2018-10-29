#ifndef CARIBOU_TOPOLOGY_GRID_H
#define CARIBOU_TOPOLOGY_GRID_H

#include <caribou/algebra/Vector.h>

namespace caribou
{
namespace topology
{

struct Grid
{
    using Vector = algebra::Vector<3, float>;

    size_t Hx; ///< Cell width
    size_t Hy; ///< Cell height
    size_t Hz; ///< Cell depth

    VectorType anchor; ///< Anchor point position. This is usually the top-left node of the grid.
};

} // namespace topology

} // namespace caribou

#endif //CARIBOU_TOPOLOGY_GRID_H
