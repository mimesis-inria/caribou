#ifndef SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_H
#define SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <Caribou/Algebra/Vector.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
#include <Caribou/Topology/Grid/Grid.h>
#include <Caribou/config.h>

#include <memory>
#include <exception>
#include <bitset>
#include <functional>

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

template <typename DataTypes>
class FictitiousGrid : public virtual BaseObject
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(FictitiousGrid, DataTypes), BaseObject);

    static constexpr unsigned char Dimension = DataTypes::spatial_dimensions;

    // Caribou data aliases
    using Index = std::size_t;
    using VecFloat = caribou::algebra::Vector<Dimension, FLOATING_POINT_TYPE>;
    using VecInt   = caribou::algebra::Vector<Dimension, Index>;
    using Int   = typename VecInt::ValueType;
    using Float = typename VecFloat::ValueType;

    // Sofa data aliases
    using SofaFloat = typename DataTypes::Real;
    using SofaVecInt = sofa::defaulttype::Vec<Dimension, Int>;
    using SofaVecFloat = sofa::defaulttype::Vec<Dimension, SofaFloat>;
    using Coord = typename DataTypes::Coord;
    using SofaVecCoord = sofa::helper::vector<Coord>;
    using ElementId = sofa::core::topology::Topology::index_type;
    using VecElementId = sofa::helper::vector<ElementId>;
    using Triange = sofa::core::topology::BaseMeshTopology::Triangle;
    using Edge = sofa::core::topology::BaseMeshTopology::Edge;

    // Grid data aliases
    using GridType = caribou::topology::Grid<Dimension>;
    using NodeIndex = typename GridType::NodeIndex;
    using CellIndex = typename GridType::CellIndex;
    using Dimensions = typename GridType::Dimensions;
    using Subdivisions = typename GridType::Subdivisions;
    using LocalCoordinates = typename GridType::LocalCoordinates;
    using WorldCoordinates = typename GridType::WorldCoordinates;
    using GridCoordinates = typename GridType::GridCoordinates;
    using CellSet = typename GridType::CellSet;
    using CellElement = typename GridType::Element;

    // FictitiousGrid aliases
    enum class Type {
        Undefined = -1,
        Inside = 0,
        Outside = 1,
        Boundary = 2
    };

    using f_implicit_test_callback_t = std::function<float(const WorldCoordinates &)>;

    template <typename ObjectType>
    using Link = SingleLink<FictitiousGrid<DataTypes>, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // Public functions
    FictitiousGrid();
    void init() override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;
    void create_grid();
    void compute_cell_types_from_implicit_surface();
    void compute_cell_types_from_explicit_surface();
    std::pair<Coord, Coord> compute_bbox_from(const SofaVecCoord & positions);

    /**
     * Set the implicit test callback function.
     * @param callback This should point to a function that takes one world position as argument and return 0 if the
     * given position is directly on the surface, < 0 if it is inside the surface, > 1 otherwise.
     *
     * float implicit_test(const WorldCoordinates & query_position);
     */
    inline void
    set_implicit_test_function(const f_implicit_test_callback_t & callback)
    {
        p_implicit_test_callback = callback;
    }

    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override
    {
        if( !onlyVisible )
            return;

        this->f_bbox.setValue(params,sofa::defaulttype::TBoundingBox<Float>(d_min.getValue().array(),d_max.getValue().array()));
    }

private:
    Data<SofaVecInt> d_n;
    Data<SofaVecFloat> d_min;
    Data<SofaVecFloat> d_max;
    Data<unsigned char> d_number_of_subdivision;
    Data<bool> d_use_implicit_surface;

    Data< SofaVecCoord > d_surface_positions;
    Data<sofa::helper::vector<Edge> > d_surface_edges; ///< List of edges (ex: [e1p1 e1p2 e2p1 e2p2 ...]).
    Data<sofa::helper::vector<Triange> > d_surface_triangles; ///< List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...]).

    std::unique_ptr<GridType> p_grid;
    f_implicit_test_callback_t p_implicit_test_callback;

    std::vector<Type> p_node_types; ///< Types of the complete regular grid's nodes
    std::vector<Type> p_cells_types; ///< Types of the complete regular grid's cells
};

template<> void FictitiousGrid<Vec2Types>::create_grid ();
template<> void FictitiousGrid<Vec3Types>::create_grid ();

template<> void FictitiousGrid<Vec2Types>::compute_cell_types_from_explicit_surface ();
template<> void FictitiousGrid<Vec3Types>::compute_cell_types_from_explicit_surface ();

extern template class FictitiousGrid<Vec2Types>;
extern template class FictitiousGrid<Vec3Types>;

} // namespace SofaCaribou::GraphComponents::topology

#endif //SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_H
