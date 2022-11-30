#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/IsoSurface.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
#include <sofa/defaulttype/Vec.h>
#endif
DISABLE_ALL_WARNINGS_END

#include <Caribou/Geometry/RectangularQuad.h>
#include <Caribou/Geometry/RectangularHexahedron.h>
#include <Caribou/Topology/Grid/Grid.h>

#include <memory>
#include <exception>
#include <bitset>
#include <functional>
#include <sstream>

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa { using Index = unsigned int; }
#endif

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210600)
namespace sofa::type {
template <typename T> using vector = ::sofa::helper::vector<T>;
using Vector3 = ::sofa::defaulttype::Vector3;
using RGBAColor = ::sofa::helper::types::RGBAColor;
template <typename Real> using TBoundingBox = ::sofa::defaulttype::TBoundingBox<Real> ;
template <std::size_t N, typename Real>
using Vec = sofa::defaulttype::Vec<N, Real>;

}
#endif

namespace SofaCaribou::topology {

using namespace sofa::core::objectmodel;

template <typename DataTypes>
class FictitiousGrid : public virtual BaseObject
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(FictitiousGrid, DataTypes), BaseObject);

    static constexpr unsigned char Dimension = DataTypes::spatial_dimensions;

    // --------------------
    // Caribou data aliases
    // --------------------
    using Index = std::size_t;
    using Int   = INTEGER_TYPE;
    using Float = FLOATING_POINT_TYPE;

    // -----------------
    // Sofa data aliases
    // -----------------
    using SofaFloat = typename DataTypes::Real;
    using SofaVecInt = sofa::type::Vec<Dimension, UNSIGNED_INTEGER_TYPE>;
    using SofaVecFloat = sofa::type::Vec<Dimension, SofaFloat>;
    using Coord = typename DataTypes::Coord;
    using SofaVecCoord = sofa::type::vector<Coord>;
    using ElementId = sofa::Index;
    using VecElementId = sofa::type::vector<ElementId>;
    using SofaHexahedron = sofa::core::topology::BaseMeshTopology::Hexahedron;
    using SofaQuad = sofa::core::topology::BaseMeshTopology::Quad;
    using SofaTriangle = sofa::core::topology::BaseMeshTopology::Triangle;
    using SofaEdge = sofa::core::topology::BaseMeshTopology::Edge;

    // -----------------
    // Grid data aliases
    // -----------------
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

    // -----------------
    // Structures
    // -----------------
    enum class Type : INTEGER_TYPE {
        Undefined = std::numeric_limits<INTEGER_TYPE>::lowest(),
        Inside =   (unsigned) 1 << (unsigned) 0,
        Outside =  (unsigned) 1 << (unsigned) 1,
        Boundary = (unsigned) 1 << (unsigned) 2
    };

    ///< The Cell data represent data attached to a leaf-cell. Only leaf-cells have this data structure assigned.
    struct CellData {
        CellData(const Type & t, const Float& w, const int & r, const bool &b)
        : type(t), weight(w), region_id(r), boundary_of_region(b) {}
        Type type = Type::Undefined;
        Float weight = 0.; // 1 for full cell, 1/4 (resp. 1/8) for the first level of subdivision in 2D (resp. 3D), etc.
        int region_id = -1;
        bool boundary_of_region = false; // True if the leaf-cell makes the boundary between its region and another one
    };

    ///< The Cell structure contains the quadtree (resp. octree) data of a given cell or subcell.
    struct Cell {
        Cell * parent = nullptr;
        CellIndex index = 0; // Index relative to the parent cell
        std::unique_ptr<CellData> data; // Data is only stored on leaf cells
        std::unique_ptr<std::array<Cell,(unsigned) 1 << Dimension>> childs;

        inline bool is_leaf() const {return childs.get() == nullptr;}
    };

    ///< A region is a cluster of cells sharing the same type and surrounded by either a boundary region or the outside
    ///< of the grid
    struct Region {
        Type type = Type::Undefined;
        std::vector<Cell*> cells; // Leaf cells filling the region
    };

    // -------
    // Aliases
    // -------

    template <typename ObjectType>
    using Link = SingleLink<FictitiousGrid<DataTypes>, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // ----------------
    // Public functions
    // ----------------
    
    FictitiousGrid();

    /** Initialization of the grid. This must be called before anything else. */
    
    virtual void
    create_grid();

    /** Get the number of sparse cells in the grid */
    inline UNSIGNED_INTEGER_TYPE
    number_of_cells() const {
        return p_cell_index_in_grid.size();
    }

    /** Get the number of sparse nodes in the grid */
    inline UNSIGNED_INTEGER_TYPE
    number_of_nodes() const {
        return p_node_index_in_grid.size();
    }

    /** Get the number of subdivisions in the grid */
    inline UNSIGNED_INTEGER_TYPE
    number_of_subdivisions() const {
        return d_number_of_subdivision.getValue();
    }

    /**
     * Get neighbors cells around a given cell. A cell is neighbor to another one if they both have a face in common,
     * or if a face contains one of the face of the other. Neighbors outside of the surface boundary are excluded.
     */
    
    std::vector<Cell *>
    get_neighbors(const Cell * cell) const;

    /**
     * Get cells neighbors in a given axis (x=0, y=1, z=2) and direction (-1, 1) of a given cell. For example,
     * to get the neighbors cells that are on the top of the given cell (in the y direction), axis would be 1 and
     * direction would be 1. Neighbors outside of the surface boundary are excluded.
     */
    
    std::vector<Cell *>
    get_neighbors(const Cell * cell, UNSIGNED_INTEGER_TYPE axis, INTEGER_TYPE direction) const;

    /**
     * Get the list of gauss nodes coordinates and their respective weight inside a cell. Here, all the gauss nodes of
     * the leafs cells that are within (or onto) the boundary are given. The coordinates are given with respect of the
     * local frame of the cell (local coordinates).
     *
     * * @param sparse_cell_index The index of the cell in the sparse grid
     */
    
    std::vector<std::pair<LocalCoordinates, FLOATING_POINT_TYPE>>
    get_gauss_nodes_of_cell(const CellIndex & sparse_cell_index) const;

    /**
     * Similar to `get_gauss_nodes_of_cell(const CellIndex & index)`, but here only the gauss nodes of inner cells up to
     * the subdivision level given are returned. Leafs cells bellow the given level are only used to compute the weight
     * of a gauss node.
     *
     * For example, if the grid's subdivision level is 3, calling this function with level = 0 will give the standard
     * 4 gauss nodes in 2D (8 gauss nodes in 3D), but where each gauss nodes will use their underlying quad tree
     * (resp. octree in 3D) to compute their weight.
     *
     * @param sparse_cell_index The index of the cell in the sparse grid
     */
    
    std::vector<std::pair<LocalCoordinates, FLOATING_POINT_TYPE>>
    get_gauss_nodes_of_cell(const CellIndex & sparse_cell_index, const UNSIGNED_INTEGER_TYPE level) const;

    /**
     * Get the element of a cell from its index in the sparse grid.
     */
    inline CellElement
    get_cell_element(const CellIndex & sparse_cell_index) const {
        const auto cell_index = p_cell_index_in_grid[sparse_cell_index];
        return std::move(p_grid->cell_at(cell_index));
    }

    /**
     * Get the node indices of a cell from its index in the sparse grid.
     */
    inline const SofaHexahedron &
    get_node_indices_of(const CellIndex & sparse_cell_index) const {
        return d_hexahedrons.getValue().at(sparse_cell_index);
    }

    /**
     * Get the type (inside, outside, boundary or undefined) of a given point in space.
     */
    
    inline Type
    get_type_at(const WorldCoordinates & p) const;

    /**
     * Compute the distribution of volume ratios of the top level cells of the grid.
     *
     * The volume ratio is the ratio of actual volume of a cell over the total volume of the cell.
     * Hence, the ratio of a cell outside the boundaries is 0, the ratio of a cell inside is 1,
     * and the ratio of boundary cells are between 0 and 1.
     *
     * @param number_of_decimals Round the volume ratio to the given
     *        number of decimals. For example, setting this value to 2  will
     *        generate a distribution of maximum 100 entries (0.00, 0.01, 0.02, ..., 0.99, 1.00).
     *
     *        Setting a value at zero deactivate the rounding of volume ratio.
     *        Default is 0 which means no rounding.
     *
     * @return A sorted map where the keys are the percentage of volume inside
     *         the cell, and the value is a vector containing the ids of all
     *         cells having this volume percentage.
     */
    
    inline std::map<FLOATING_POINT_TYPE, std::vector<CellIndex>>
    cell_volume_ratio_distribution(UNSIGNED_INTEGER_TYPE number_of_decimals=0) const;

    // ---------------------
    // SOFA METHOD OVERRIDES
    // ---------------------
    
    void init() override;

    
    void draw(const sofa::core::visual::VisualParams* vparams) override;

    void computeBBox(const sofa::core::ExecParams*, bool onlyVisible) override {
        if( !onlyVisible )
            return;

        if (Dimension == 2) {
            const Float min[3] = {
                d_min.getValue()[0], d_min.getValue()[1], -1
            };
            const Float max[3] = {
                d_max.getValue()[0], d_max.getValue()[1], +1
            };
            this->f_bbox.setValue(sofa::type::TBoundingBox<Float>(min, max));
        } else {
            this->f_bbox.setValue(sofa::type::TBoundingBox<Float>(
                d_min.getValue().array(),d_max.getValue().array()));
        }
    }

    static std::string GetCustomTemplateName() {
        return templateName();
    }

    static std::string templateName(const FictitiousGrid<DataTypes>* = nullptr) {
        return DataTypes::Name();
    }

private:
    virtual void tag_intersected_cells_from_implicit_surface();
    virtual void tag_intersected_cells();
    virtual void tag_outside_cells();
    virtual void tag_inside_cells();
    virtual void subdivide_intersected_cells();
    virtual void create_regions_from_same_type_cells();
    virtual void create_sparse_grid();
    virtual void populate_drawing_vectors();
    virtual void validate_grid();

    std::array<CellElement, (unsigned) 1 << Dimension> get_subcells_elements(const CellElement & e) const;
    std::vector<Cell *> get_leaf_cells(const Cell & c) const {return std::move(get_leaf_cells(&c));}
    std::vector<Cell *> get_leaf_cells(const Cell * c) const;
    inline FLOATING_POINT_TYPE get_cell_weight(const Cell & c) const;

private:
    // ------------------
    // Input data members
    // ------------------
    Data<SofaVecInt> d_n;
    Data<SofaVecFloat> d_min;
    Data<SofaVecFloat> d_max;
    Data<UNSIGNED_INTEGER_TYPE> d_number_of_subdivision;
    Data<Float> d_volume_threshold;
    Link<IsoSurface<DataTypes>> d_iso_surface;
    Data<bool> d_draw_boundary_cells;
    Data<bool> d_draw_outside_cells;
    Data<bool> d_draw_inside_cells;
    Data<double> d_draw_scale;

    Data< SofaVecCoord > d_surface_positions;

    ///< List of edges (ex: [e1p1 e1p2 e2p1 e2p2 ...]).
    Data<sofa::type::vector<SofaEdge> > d_surface_edges;

    ///< List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...]).
    Data<sofa::type::vector<SofaTriangle> > d_surface_triangles;

    // -------------------
    // Output data members
    // -------------------
    ///< Position vector of nodes contained in the sparse grid
    Data< SofaVecCoord > d_positions;

    ///< List of quads contained in the sparse grid (ex: [q1p1 q1p2 q1p3 q1p4 q2p1 ... qnp3 qnp4]).
    Data < sofa::type::vector<SofaQuad> > d_quads;

    ///< List of hexahedrons contained in the sparse grid (ex: [h1p1 h1p2 h1p3 h1p4 h1p5 ... hnp6 hnp7]).
    Data < sofa::type::vector<SofaHexahedron> > d_hexahedrons;

    // ---------------
    // Private members
    // ---------------
    ///< The underground grid object. This object do not store any values beside the dimensions and size of the grid.
    ///< Most of the grid algorithms are defined there.
    std::unique_ptr<GridType> p_grid;

    ///< Types of the complete regular grid's cells
    std::vector<Type> p_cells_types;

    ///< List of boundary elements that intersect a given cell.
    std::vector<std::vector<Index>> p_triangles_of_cell;

    ///< Quadtree (resp. Octree) representation of the 2D (resp 3D) cell.
    std::vector<Cell> p_cells;

    ///< Distinct regions of cells.
    std::vector<Region> p_regions;

    ///< Contains the index of a node in the sparse grid from its index in the full grid, or -1 if the node isn't
    ///< present in the sparse grid.
    std::vector<INTEGER_TYPE> p_node_index_in_sparse_grid;

    ///< Contains the index of a node in the full grid from its index in the sparse grid.
    std::vector<UNSIGNED_INTEGER_TYPE> p_node_index_in_grid;

    ///< Contains the index of a cell in the sparse grid from its index in the full grid, or -1 if the cell isn't
    ///< present in the sparse grid.
    std::vector<INTEGER_TYPE> p_cell_index_in_sparse_grid;

    ///< Contains the index of a cell in the full grid from its index in the sparse grid.
    std::vector<UNSIGNED_INTEGER_TYPE> p_cell_index_in_grid;

    ///< Contains the grid's nodes to be draw
    std::vector<sofa::type::Vector3> p_drawing_nodes_vector;

    ///< Contains the grid's edges to be draw
    std::vector<sofa::type::Vector3> p_drawing_edges_vector;

    ///< Contains the edges of subdivided cells for each region to be draw
    std::vector<std::vector<sofa::type::Vector3>> p_drawing_subdivided_edges_vector;

    ///< Contains the cells for each region to be draw
    std::vector<std::vector<sofa::type::Vector3>> p_drawing_cells_vector;

    // ----------------------
    // Private static members
    // ----------------------
    ///< Contains the coordinates of each subcells of a quad (resp hexa) in 2D resp(3D)
    static const GridCoordinates subcell_coordinates[(unsigned) 1 << Dimension];

};

#ifndef WIN32
template<> const FictitiousGrid<sofa::defaulttype::Vec2Types>::GridCoordinates FictitiousGrid<sofa::defaulttype::Vec2Types>::subcell_coordinates[4];
template<> const FictitiousGrid<sofa::defaulttype::Vec3Types>::GridCoordinates FictitiousGrid<sofa::defaulttype::Vec3Types>::subcell_coordinates[8];
#endif

template<> void FictitiousGrid<sofa::defaulttype::Vec2Types>::tag_intersected_cells ();
template<> void FictitiousGrid<sofa::defaulttype::Vec3Types>::tag_intersected_cells ();

template<> void FictitiousGrid<sofa::defaulttype::Vec2Types>::subdivide_intersected_cells ();
template<> void FictitiousGrid<sofa::defaulttype::Vec3Types>::subdivide_intersected_cells ();

extern template class FictitiousGrid<sofa::defaulttype::Vec2Types>;
extern template class FictitiousGrid<sofa::defaulttype::Vec3Types>;

} // namespace SofaCaribou::topology
