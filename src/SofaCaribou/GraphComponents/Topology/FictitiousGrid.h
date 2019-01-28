#ifndef SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_H
#define SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>

#include <Caribou/Algebra/Vector.h>
#include <Caribou/Topology/Engine/Grid/Grid.h>
#include <Caribou/config.h>

#include <memory>

// Forward declarations
namespace caribou::topology::engine {
template <size_t Dim>
struct Grid;
}

namespace sofa::helper::io {
struct Mesh;
}

namespace SofaCaribou::GraphComponents::topology {

using namespace sofa::core::objectmodel;

using sofa::defaulttype::Vec2Types;
using sofa::defaulttype::Vec3Types;

template <typename VecType>
class FictitiousGrid : public virtual BaseObject
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(FictitiousGrid, VecType), BaseObject);

    static constexpr char Dimension = VecType::spatial_dimensions;
    using GridType = caribou::topology::engine::Grid<Dimension>;

    // Grid data aliases
    using NodeIndex = typename GridType::NodeIndex;
    using CellIndex = typename GridType::CellIndex;
    using Dimensions = typename GridType::Dimensions;
    using Subdivisions = typename GridType::Subdivisions;
    using LocalCoordinates = typename GridType::LocalCoordinates;
    using WorldCoordinates = typename GridType::WorldCoordinates;
    using GridCoordinates = typename GridType::GridCoordinates;
    using CellSet = typename GridType::CellSet;

    // Caribou data aliases
    using Index = size_t;
    using VecFloat = caribou::algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using VecInt   = caribou::algebra::Vector<Dimension, Index>;
    using Int   = typename VecInt::ValueType;
    using Float = typename VecFloat::ValueType;

    // Sofa data aliases
    using SofaFloat = typename VecType::Real;
    using SofaVecInt = sofa::defaulttype::Vec<Dimension, Int>;
    using SofaVecFloat = sofa::defaulttype::Vec<Dimension, SofaFloat>;
    using Coord = typename VecType::Coord;

    // Immersed domain data aliases
    using DomainIndex = Index;

    /**
     * An immersed domain representation
     */
    struct ImmersedDomain {
        std::vector<Coord> boundary_positions;
        std::vector<Index> boundary_element_indices;
    };

    /**
     * One Cell representation
     */
    struct Cell {
        ///< Each cell can be intersected by one, or even multiple immersed domains
        std::vector<DomainIndex> immersed_domains;
    };

    FictitiousGrid();
    void init() override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;

private:
    Data<SofaVecInt> d_n;
    Data<SofaVecFloat> d_min;
    Data<SofaVecFloat> d_max;
    Data<unsigned char> d_number_of_subdivision;

    std::unique_ptr<GridType> p_grid;
    std::vector<Cell> p_cells;

};

extern template class FictitiousGrid<Vec2Types>;
extern template class FictitiousGrid<Vec3Types>;

} // namespace SofaCaribou::GraphComponents::topology

#endif //SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_H
