#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Vector.h>
#include <Caribou/Topology/Engine/Grid/Internal/BaseCell.h>
#include <memory>

namespace caribou::topology::engine {

/**
 * A Cell represents a segment in 1D, a rectangle quad in 2D or a 3D rectangle hexahedron in 3D.
 *
 * @tparam Dim The dimension of the cell (1D, 2D or 3D).
 */
template <unsigned char Dim>
struct Cell : public internal::BaseCell<Dim, Cell<Dim>>
{
    static constexpr char Dimension = Dim;
    static_assert(Dimension == 1 or Dimension == 2 or Dimension == 3, "A cell must be in one, two or three dimension.");
};

/** One dimensional cell (segment) **/
template <>
struct Cell<1> : public internal::BaseCell<1, Cell<1>>
{
    using Base = internal::BaseCell<1, Cell<1>>;

    static constexpr char Dimension = 1;
    static constexpr size_t NumberOfNodes = (unsigned char) (1 << Dimension);
};

/** Two dimensional cell (square) **/
template <>
struct Cell<2> : public internal::BaseCell<2, Cell<2>>
{
    using Base = internal::BaseCell<2, Cell<2>>;

    static constexpr char Dimension = 2;
    static constexpr size_t NumberOfNodes = (unsigned char) (1 << Dimension);
};

/** Three dimensional cell (cube) **/
template <>
struct Cell<3> : public internal::BaseCell<3, Cell<3>>
{
    using Base = internal::BaseCell<3, Cell<3>>;

    static constexpr char Dimension = 3;
    static constexpr size_t NumberOfNodes = (unsigned char) (1 << Dimension);
};

} // namespace caribou::topology::engine

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_CELL_H
