#include <Caribou/Topology/Grid.h>
#include <exception>

namespace caribou
{

using namespace geometry;

namespace topology
{

using VecFloat = Grid::VecFloat;
using VecInt = Grid::VecInt;

Grid::Grid(VecFloat anchor, VecInt subdivisions, VecFloat dimensions)
        : parent(nullptr), anchor(anchor), dimensions(dimensions), cells()
{
    subdivide(subdivisions);
};

Grid::Grid(Grid* parent, VecFloat anchor, VecFloat dimensions)
        : parent(parent), anchor(anchor), dimensions(dimensions), cells()
{
};

void Grid::subdivide(VecInt subdivisions)
{
    if (!is_a_leaf()) {
        throw std::logic_error("Trying to subdivide an already subdivided grid.");
    }

    const auto & nx = subdivisions[0];
    const auto & ny = subdivisions[1];
    const auto & nz = subdivisions[2];

    const auto h = cell_size(subdivisions);

    cells.resize(nz);
    for (size_t k = 0; k < nz; ++k) {
        cells[k].resize(ny);
        for (size_t j = 0; j < ny; ++j) {
            cells[k][j].resize(nx);
            for (size_t i = 0; i < nx; ++i) {
                cells[k][j][i] =
                    std::make_unique<Grid>(
                            this,
                            cell_anchor({i, j, k}), /* Anchor position in local coordinate */
                            h /* Size of the cell to subdivide */
                    );
            }
        }
    }
};

VecInt Grid::subdivisions() const
{
    if (is_a_leaf())
        return {0, 0, 0};

    const auto & nz = cells.size();
    const auto & ny = cells[0].size();
    const auto & nx = cells[0][0].size();

    return {nx, ny, nz};
};

VecFloat Grid::cell_size() const
{
    return cell_size(subdivisions());
};


VecFloat Grid::cell_size(VecInt subdivisions) const
{
    const auto & n = subdivisions;
    const auto & nx = n[0];
    const auto & ny = n[1];
    const auto & nz = n[2];

    const auto & s = size();
    const auto & sx = s[0];
    const auto & sy = s[1];
    const auto & sz = s[2];

    const auto hx = (nx==0) ? sx : sx/nx;
    const auto hy = (ny==0) ? sy : sy/ny;
    const auto hz = (nz==0) ? sz : sz/nz;

    return {hx, hy, hz};
};

VecFloat Grid::cell_anchor(const VecInt & grid_coordinates) const
{
    const VecInt::ValueType & i = grid_coordinates[0];
    const VecInt::ValueType & j = grid_coordinates[1];
    const VecInt::ValueType & k = grid_coordinates[2];

    const auto & n = subdivisions();
    const auto & nx = n[0];
    const auto & ny = n[1];
    const auto & nz = n[2];

    if (i > nx-1 || j > ny-1 || k > nz) {
        throw std::out_of_range(
            "Trying to access a cell at an invalid grid coordinate (" + std::to_string(i) + ", " + std::to_string(j) + "," + std::to_string(j) + ")"
        );
    }

    return {
            (2. * (float)i / (float) nx) - 1., // xi
            (2. * (float)j / (float) ny) - 1., // eta
            (2. * (float)k / (float) nz) - 1., // zeta
    };
};

Grid & Grid::cell(const VecInt & grid_coordinates)
{
    const VecInt::ValueType & i = grid_coordinates[0];
    const VecInt::ValueType & j = grid_coordinates[1];
    const VecInt::ValueType & k = grid_coordinates[2];

    const auto & n = subdivisions();
    const auto & nx = n[0];
    const auto & ny = n[1];
    const auto & nz = n[2];

    if (i > nx-1 || j > ny-1 || k > nz) {
        throw std::out_of_range(
                "Trying to access a cell at an invalid grid coordinate (" + std::to_string(i) + ", " + std::to_string(j) + "," + std::to_string(j) + ")"
        );
    }

    return *cells[k][j][i];
};

VecFloat Grid::position(const VecFloat & local_coordinate) const
{
    const auto & s = size();
    auto sx = s[0];
    auto sy = s[1];
    auto sz = s[2];

    if (parent) {
        const auto & h = parent->size();
        sx = sx / h[0] * 2;
        sy = sy / h[1] * 2;
        sz = sz / h[2] * 2;
    }

    const VecFloat coordinates = {
            anchor[0] + 1/2. * (1 + local_coordinate[0])*sx,
            anchor[1] + 1/2. * (1 + local_coordinate[1])*sy,
            anchor[2] + 1/2. * (1 + local_coordinate[2])*sz
    };

    if (parent) {
        return parent->position(coordinates);
    } else {
        return coordinates;
    }
}

//std::vector<VecFloat>
//Grid::positions(bool include_cells) const
//{
//
//};

} // namespace topology

} // namespace caribou