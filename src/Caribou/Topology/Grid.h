#ifndef CARIBOU_TOPOLOGY_GRID_H
#define CARIBOU_TOPOLOGY_GRID_H

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

struct Grid
{
    using VecFloat = algebra::Vector<3, FLOATING_POINT_TYPE>;
    using VecInt = algebra::Vector<3, unsigned int>;

    /** Default constructor is not permitted **/
    Grid() = delete;

    /**
     * Constructor of the top level grid (the parent of all subdivisions)
     * @param anchor The anchor point position vector (x, y, z). This is the position of the node id #0 on a hexahedron (see caribou::geometry::LinearHexahedron).
     * @param subdivisions Vector of 3 integer (nx, ny, nz) which specify the number of subcells in the x, y and z directions respectively
     * @param dimensions Vector of 3 float (sx, sy, sz) which specify the dimension of the grid from the anchor point in the x, y and z directions respectively
     * @param parent If null, the grid will be initialized as the top level grid. Else, the grid is a subcell of the Grid parent.
     */
    Grid(VecFloat anchor, VecInt subdivisions, VecFloat dimensions);

    /**
     * Constructor of a sub-cell grid.
     * @param parent Parent grid of this new grid.
     * @param anchor The anchor point position vector (x, y, z). This is the position of the node id #0 on a hexahedron (see caribou::geometry::LinearHexahedron).
     * @param dimensions Vector of 3 float (sx, sy, sz) which specify the dimension of the grid from the anchor point in the x, y and z directions respectively
     */
    Grid(Grid* parent, VecFloat anchor, VecFloat dimensions);

    /**
     * Subdivide the Grid into nx, ny and nz subcells
     * @param subdivisions Specify the number (nx, ny, nz) of subcells
     * @throws std::logic_error When the grid is already subdivided into a number of cells.
     */
    void subdivide(VecInt subdivisions = {2,2,2});

    /**
     * Compute the local coordinates of the anchor point of the cell located at grid_index (i, j, k).
     *
     * The anchor point is located on the node id #0 (see caribou::geometry::LinearHexahedron::node) of the cell.
     *
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     *
     * @return  The local coordinates (xi, eta, zeta) of the cell's anchor point relative to the grid anchor point.
     *
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    VecFloat cell_anchor(const VecInt & grid_coordinates) const;

    /** Get this grid dimensions (sx, sy, sz). **/
    inline VecFloat size() const {return dimensions;}

    /** Get dimensions (hx, hy, hz) of a cell in this grid. **/
    VecFloat cell_size() const;

    /** Compute the cell size if we divide the current grid into (nx, ny, nz) subdivisions. **/
    VecFloat cell_size(VecInt subdivisions) const;

    /** Get the number of cell subdivisions (nx, ny, nz) of this grid. **/
    VecInt subdivisions() const;

    /**
     * Get the cell located at grid_index (i, j, k).
     * @param grid_coordinates Cell location provided in terms of grid coordinates (i, j, k).
     * @throws std::out_of_range when the grid coordinates are outside of this grid subdivisions (nx, ny, nz).
     */
    Grid & cell(const VecInt & grid_coordinates);

    /** Alias to caribou::topology::cell(const VecInt & grid_coordinates) **/
    inline Grid & cell(const VecInt::ValueType i, const VecInt::ValueType j, const VecInt::ValueType k) {return cell({i, j, k});};

    /** True if the grid is a leaf (it contains no sub-cells) **/
    inline bool is_a_leaf() const {return cells.empty();};

    /** Compute the global position (x,y,z) of a point in this grid - in terms of local coordinate (xi, eta, zeta).  **/
    VecFloat position(const VecFloat & local_coordinate) const;

    /**
     * Get the positions vector of this grid.
     * @param include_subcells If true, the positions of the subcells (recursively) will be added to the vector. Note that this won't produce duplicated nodes.
     */
//    std::vector<VecFloat> positions(bool include_subcells = false) const;

protected:
    Grid* parent; ///< The parent grid of this grid. It should be null for the top-level grid.
    VecFloat anchor; ///< The anchor point position. This is the position of the node id #0 on a hexahedron.
    VecFloat dimensions; ///< The anchor point position. This is the position of the node id #0 on a hexahedron.
    std::vector<std::vector<std::vector<std::unique_ptr<Grid>>>> cells; ///< This grid can be subdivided into a set of sub-cells where a sub-cell can be a leaf (no further subdivisions) or a grid.
//    size_t level; ///< The level of this grid in the subdivision graph. It starts with 0 on the top-most grid. If this grid's level is n, its subcells levels are (n+1).
};

} // namespace topology

} // namespace caribou

#endif //CARIBOU_TOPOLOGY_GRID_H
