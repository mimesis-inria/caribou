#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H

#include <vector>
#include <memory>

#include <Caribou/config.h>

namespace caribou
{

namespace topology
{

namespace engine
{

/**
 * A Grid is a rectangular 2D quad (resp. 3D hexahedron) that contain multiple Cell entities aligned in the x, y (and z in 3D) axis.
 * @tparam Dim The dimension (2D or 3D) of the grid.
 * @tparam TCell The type of cell this grid will use.
 */
template <unsigned char Dim, class TCell>
struct BaseGrid
{
    using CellType = TCell;
    using VecFloat = typename CellType::VecFloat;
    using Index = typename CellType::Index;
    using VecInt = typename CellType::VecInt;

    static constexpr size_t Dimension = Dim;
    static constexpr size_t NumberOfNodes = CellType::NumberOfNodes;

    static_assert(Dimension == CellType::Dimension, "The dimension of the grid must match the dimension of its cell type.");

    /** Default constructor is not permitted **/
    BaseGrid() = delete;

    /**
     * Constructor of a grid
     * @param anchor The anchor point position vector (x, y, z). It should be positioned at the center of the grid.
     * @param subdivisions Vector of 3 integer (nx, ny, nz) which specify the number of sub-cells in the x, y and z directions respectively
     * @param dimensions Vector of 3 float (sx, sy, sz) which specify the dimension of the grid from the anchor point in the x, y and z directions respectively
     * @param parent If null, the grid will be initialized as the top level grid. Else, the grid is a subcell of the Grid parent.
     */
    BaseGrid(VecFloat /*anchor*/, VecInt /*subdivisions*/, VecFloat /*dimensions*/) {}

    /** Get the number of cell subdivisions (nx, ny, nz) of this grid. **/
    inline VecInt
    number_of_subdivision() const {return nSubdivisions;};

    /** Get this grid dimensions (sx, sy, sz). **/
    inline VecFloat
    size() const { return dimensions;};

    /** Get dimensions (hx, hy, hz) of a cell in this grid. **/
    inline VecFloat
    cell_size() const
    {
        const auto & N = number_of_subdivision();
        const auto & S = size();
        const auto   H = S.direct_division(N);

        return H;
    };

    /**
     * Get the cell index at grid location (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    Index cell_index(const VecInt & /*grid_coordinates*/) const {return (Index) -1;}

    /** Get the grid location (i, j, k) at index cell_index */
    VecInt grid_coordinates(const Index & /*cell_index*/) const {return {-1, -1, -1};};

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    inline CellType &
    get(const VecInt & grid_coordinates)
    {
        return *(cells[cell_index(grid_coordinates)]);
    }

    /**
     * Get the indices of the nodes forming the cell located at grid_index (i, j, k). The indices of the nodes will be return
     * following caribou::geometry::LinearHexahedron::node node numbering.
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    std::array<Index, NumberOfNodes> nodes(const VecInt & grid_coordinates) const;

    /** Get the position of the node node_id **/
    VecFloat position(const Index & node_id) const;

protected:
    ///< The anchor point position. It should be positioned at the center of the grid.
    VecFloat anchor;

    ///< Number of sub-cells in the x, y and z directions respectively.
    const VecInt nSubdivisions;

    ///<  Dimension of the grid from the anchor point in the x, y and z directions respectively.
    const VecFloat dimensions;

    ///< This grid can be subdivided into a set of sub-cells where a sub-cell can be a leaf (no further subdivisions) or a grid.
    std::vector<std::unique_ptr<CellType>> cells;
};

template <unsigned char Dim, class TCell>
struct Grid : public BaseGrid<Dim, TCell>
{
};

/** 2D partial specialization of caribou::topology::engine::BaseGrid */
template <class TCell>
struct Grid<2, TCell> : public BaseGrid<2, TCell>
{
    using Base = BaseGrid<2, TCell>;
    using typename Base::VecFloat;
    using typename Base::VecInt;
    using typename Base::Index;
    using Base::NumberOfNodes;

    /** see caribou::topology::engine::BaseGrid::BaseGrid */
    Grid(VecFloat anchor, VecInt subdivisions, VecFloat dimensions);

    /** see caribou::topology::engine::BaseGrid::cell_index */
    Index cell_index(const VecInt & grid_coordinates) const;

    /** see caribou::topology::engine::BaseGrid::grid_coordinates */
    VecInt grid_coordinates(const Index & cell_index) const;

    /** see caribou::topology::engine::BaseGrid::nodes */
    std::array<Index, NumberOfNodes> nodes(const VecInt & grid_coordinates) const;
};

/** 3D partial specialization of caribou::topology::engine::BaseGrid */
template <class TCell>
struct Grid<3, TCell> : public BaseGrid<3, TCell>
{
    // Importing the BaseGrid constructors
    using BaseGrid<3, TCell>::BaseGrid;
};

template <class TCell>
using Grid2D = Grid<2, TCell>;

template <class TCell>
using Grid3D = Grid<3, TCell>;

} // namespace engine

} // namespace topology

} // namespace caribou

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
