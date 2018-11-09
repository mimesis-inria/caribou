#ifndef CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
#define CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H

#include <vector>
#include <memory>

#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Hexahedron.h>

#ifdef CARIBOU_USE_DOUBLE
#define FLOATING_POINT_TYPE double
#else
#define FLOATING_POINT_TYPE float
#endif

namespace caribou
{

using namespace geometry;

namespace topology
{

namespace engine
{

template <unsigned char Dimension>
struct Grid;

/**
 * A Cell represents a 2D quad (resp. 3D hexahedron) entity in a rectangular Grid made of multiple cells. It is always
 * contained in a parent grid, and it can also be subdivided in sub-cells by adding it a grid. Otherwise, if no grid is
 * added, the cell is said to be a leaf of its parent grid.
 *
 * @tparam Dimension The dimension (2D or 3D) of the cell.
 */
template <unsigned char Dimension>
struct Cell
{
    static_assert(Dimension == 2 or Dimension == 3, "A grid cell must be in two or three dimension.");

    using VecFloat = algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using Index = size_t;
    using VecInt = algebra::Vector<Dimension, Index>;
    using GridType = Grid<Dimension>;

    static constexpr size_t NumberOfNodes = (unsigned char) (1 << Dimension);

    explicit Cell() = delete;

    Cell(Grid<Dimension>* parent, Index index) : m_grid(nullptr), m_parent(parent), m_index (index) {
        if (!parent) {
            throw std::logic_error("A cell must be contained in a parent grid.");
        }
    };

    /**
     * Subdivide the cell into nx, ny and nz sub-cells
     * @param subdivisions Specify the number (nx, ny, nz) of sub-cells
     * @throws std::logic_error When the cell is already subdivided.
     * @return The created grid containing the sub-cells
     */
    GridType * subdivide(VecInt subdivisions);

    /** This cell size (hx, hy, hz) **/
    VecFloat size() const;

    /** True if the cell is a leaf (it contains no sub-cells) **/
    inline bool
    is_a_leaf() const
    {return (m_grid == nullptr);}

    /** Get the index of the cell (relative to its parent grid). **/
    inline Index index() const {return m_index;};

    /** Get the indices of the nodes forming the cell. */
    std::array<Index, NumberOfNodes> nodes() const;

protected:
    std::unique_ptr<GridType> m_grid = nullptr;
    GridType* m_parent = nullptr;
    const Index m_index;
};

/**
 * A Grid is a rectangular 2D quad (resp. 3D hexahedron) that contain multiple Cell entities aligned in the x, y (and z in 3D) axis.
 * @tparam Dimension The dimension (2D or 3D) of the grid.
 */
template <unsigned char Dimension>
struct Grid
{
    static_assert(Dimension == 2 or Dimension == 3, "A grid must be in two or three dimension.");

    using CellType = Cell<Dimension>;
    using VecFloat = typename CellType::VecFloat;
    using Index = typename CellType::Index;
    using VecInt = typename CellType::VecInt;

    static constexpr size_t NumberOfNodes = CellType::NumberOfNodes;

    /** Default constructor is not permitted **/
    Grid() = delete;

    /**
     * Constructor of a grid
     * @param anchor The anchor point position vector (x, y, z). It should be positioned at the center of the grid.
     * @param subdivisions Vector of 3 integer (nx, ny, nz) which specify the number of sub-cells in the x, y and z directions respectively
     * @param dimensions Vector of 3 float (sx, sy, sz) which specify the dimension of the grid from the anchor point in the x, y and z directions respectively
     * @param parent If null, the grid will be initialized as the top level grid. Else, the grid is a subcell of the Grid parent.
     */
    Grid(VecFloat anchor, VecInt subdivisions, VecFloat dimensions);

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
    Index cell_index(const VecInt & grid_coordinates) const;

    /** Get the grid location (i, j, k) at index cell_index */
    VecInt grid_coordinates(const Index & cell_index) const;

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    CellType & get(const VecInt & grid_coordinates);

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


using Grid2D = Grid<2>;
using Grid3D = Grid<3>;

} // namespace engine

} // namespace topology

} // namespace caribou

#include <Caribou/Topology/Engine/Grid/Grid.inl>

#endif //CARIBOU_TOPOLOGY_ENGINE_GRID_GRID_H
