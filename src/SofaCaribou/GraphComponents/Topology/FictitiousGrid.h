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
    using Triangle = sofa::core::topology::BaseMeshTopology::Triangle;
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

    // Structures
    enum class Type : INTEGER_TYPE {
        Undefined = -1,
        Inside = 0,
        Outside = 1,
        Boundary = 2
    };

    ///< The Cell structure contains the quadtree (resp. octree) data of a given cell or subcell.
    struct CellData {
        CellData(const Type & t, const Float& w, const int & r)
        : type(t), weight(w), region_id(r) {}
        Type type = Type::Undefined;
        Float weight = 0.;
        int region_id = -1;
    };

    struct Cell {
        Cell * parent = nullptr;
        CellIndex index = 0; // Index relative to the parent cell
        std::unique_ptr<CellData> data; // Data is only stored on leaf cells
        std::unique_ptr<std::array<Cell,(unsigned) 1 << Dimension>> childs;
    };

    ///< A region is a cluster of cells sharing the same type and surrounded by either a boundary region or the outside
    ///< of the grid
    struct Region {
        Type type = Type::Undefined;
        std::vector<Cell*> cells;
    };

    // Aliases
    using f_implicit_test_callback_t = std::function<float(const WorldCoordinates &)>;

    template <typename ObjectType>
    using Link = SingleLink<FictitiousGrid<DataTypes>, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // Public functions
    FictitiousGrid();
    void init() override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;
    virtual void create_grid();
    std::vector<Cell *> get_neighbors(Cell * cell);
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

        this->f_bbox.setValue(params,sofa::defaulttype::TBoundingBox<Float>(
            d_min.getValue().array(),d_max.getValue().array()));
    }

private:
    virtual void compute_cell_types_from_implicit_surface();
    virtual void compute_cell_types_from_explicit_surface();
    virtual void subdivide_cells();
    virtual void populate_drawing_vectors();

    std::array<CellElement, (unsigned) 1 << Dimension> get_subcells(const CellElement & e) const;

private:
    // ------------
    // Data members
    // ------------
    Data<SofaVecInt> d_n;
    Data<SofaVecFloat> d_min;
    Data<SofaVecFloat> d_max;
    Data<UNSIGNED_INTEGER_TYPE> d_number_of_subdivision;
    Data<bool> d_use_implicit_surface;
    Data<bool> d_draw_boundary_cells;
    Data<bool> d_draw_outside_cells;
    Data<bool> d_draw_inside_cells;

    Data< SofaVecCoord > d_surface_positions;

    ///< List of edges (ex: [e1p1 e1p2 e2p1 e2p2 ...]).
    Data<sofa::helper::vector<Edge> > d_surface_edges;

    ///< List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...]).
    Data<sofa::helper::vector<Triangle> > d_surface_triangles;

    // ---------------
    // Private members
    // ---------------
    ///< The underground grid object. This object do not store any values beside the dimensions and size of the grid.
    ///< Most of the grid algorithms are defined there.
    std::unique_ptr<GridType> p_grid;

    ///< This is a pointer to a callback function that determines if a position is inside, outside or on the boundary.
    ///< It is used when an implicit surface definition is avaible.
    f_implicit_test_callback_t p_implicit_test_callback;

    ///< Types of the complete regular grid's nodes
    std::vector<Type> p_node_types;

    ///< Types of the complete regular grid's cells
    std::vector<Type> p_cells_types;

    ///< List of boundary elements that intersect a given cell.
    std::vector<std::vector<Index>> p_triangles_of_cell;

    ///< Quadtree (resp. Octree) representation of the 2D (resp 3D) cell.
    std::vector<Cell> p_cells;

    ///< Distinct regions of cells.
    std::vector<Region> p_regions;

    ///< Contains the grid's nodes to be draw
    std::vector<sofa::defaulttype::Vector3> p_drawing_nodes_vector;

    ///< Contains the grid's edges to be draw
    std::vector<sofa::defaulttype::Vector3> p_drawing_edges_vector;

    ///< Contains the edges of subdivided cells for each region to be draw
    std::vector<std::vector<sofa::defaulttype::Vector3>> p_drawing_subdivided_edges_vector;

    ///< Contains the cells for each region to be draw
    std::vector<std::vector<sofa::defaulttype::Vector3>> p_drawing_cells_vector;

    // ----------------------
    // Private static members
    // ----------------------
    ///< Contains the coordinates of each subcells of a quad (resp hexa) in 2D resp(3D)
    static const GridCoordinates subcell_coordinates[(unsigned) 1 << Dimension];

};

//template<> const FictitiousGrid<Vec2Types>::GridCoordinates FictitiousGrid<Vec2Types>::subcell_coordinates[4];
template<>  const FictitiousGrid<Vec3Types>::GridCoordinates FictitiousGrid<Vec3Types>::subcell_coordinates[8];

//template<> void FictitiousGrid<Vec2Types>::create_grid ();
template<> void FictitiousGrid<Vec3Types>::create_grid ();

//template<> void FictitiousGrid<Vec2Types>::compute_cell_types_from_explicit_surface ();
template<> void FictitiousGrid<Vec3Types>::compute_cell_types_from_explicit_surface ();

//template<> void FictitiousGrid<Vec2Types>::subdivide_cells ();
template<> void FictitiousGrid<Vec3Types>::subdivide_cells ();

//template<> std::array<FictitiousGrid<Vec2Types>::CellElement, (unsigned) 1 << FictitiousGrid<Vec2Types>::Dimension> FictitiousGrid<Vec2Types>::
//    get_subcells(const CellElement & e) const;

template<> std::array<FictitiousGrid<Vec3Types>::CellElement, (unsigned) 1 << FictitiousGrid<Vec3Types>::Dimension> FictitiousGrid<Vec3Types>::
get_subcells(const CellElement & e) const;

template<> void FictitiousGrid<Vec3Types>::draw (const sofa::core::visual::VisualParams* vparams);

//extern template class FictitiousGrid<Vec2Types>;
extern template class FictitiousGrid<Vec3Types>;

} // namespace SofaCaribou::GraphComponents::topology

#endif //SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_H
