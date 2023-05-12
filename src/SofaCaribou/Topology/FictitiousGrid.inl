#pragma once

#include <SofaCaribou/Topology/FictitiousGrid.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/topology/container/constant/MeshTopology.h>
#include <sofa/simulation/Node.h>
#include <sofa/core/behavior/MechanicalState.h>
DISABLE_ALL_WARNINGS_END

#include <stack>
#include <queue>
#include <iomanip>
#include <chrono>
#include <memory>
#include <tuple>
#include <bitset>

#define BEGIN_CLOCK ;std::chrono::steady_clock::time_point __time_point_begin;
#define TICK ;__time_point_begin = std::chrono::steady_clock::now();
#define TOCK (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - __time_point_begin).count())


namespace SofaCaribou::topology {

template <typename DataTypes>
FictitiousGrid<DataTypes>::FictitiousGrid()
        : d_n(initData(&d_n, SofaVecInt(),"n","Grid resolution [nx, ny, nz] where nx, ny or nz are the number of nodes in the x,y and z directions."))
        , d_min(initData(&d_min, SofaVecFloat(),"min","First corner node position of the grid's bounding box."))
        , d_max(initData(&d_max, SofaVecFloat(),"max","Second corner node position of the grid's bounding box."))
        , d_number_of_subdivision(initData(&d_number_of_subdivision,
                (UNSIGNED_INTEGER_TYPE) 0,
                "maximum_number_of_subdivision_levels",
                "Number of subdivision levels of the boundary cells (one level split the cell in 4 subcells in 2D, and 8 subcells in 3D)."))
        , d_volume_threshold(initData(&d_volume_threshold, (Float) 0.0, "volume_threshold", "Ignore every cells having a volume ratio smaller than this threshold."))
        , d_iso_surface(initLink(
                "iso_surface",
                "Use an implicit surface instead of a tessellated surface. This will be used as a level-set where an iso-value less than zero means the point is inside the boundaries."))
        , d_draw_boundary_cells(initData(&d_draw_boundary_cells,
                bool(false),
                "draw_boundary_cells",
                "Draw the cells intersected by the surface boundary."))
        , d_draw_outside_cells(initData(&d_draw_outside_cells,
                bool(false),
                "draw_outside_cells",
                "Draw the cells that are outside of the surface boundary."))
        , d_draw_inside_cells(initData(&d_draw_inside_cells,
                bool(false),
                "draw_inside_cells",
                "Draw the cells that are inside of the surface boundary."))
        , d_draw_scale(initData(&d_draw_scale,
                0.85,
                "draw_scale",
                "Scaling factor for the drawing of elements (between 0 and 1). The factor allows to shrink "
                "the element relative to its center point when drawing it."))
        , d_surface_positions(initData(&d_surface_positions, SofaVecCoord(),
                "surface_positions",
                "Position vector of nodes contained in the surface boundary of the immersed object."))
        , d_surface_edges(initData(&d_surface_edges,
                "surface_edges",
                "List of edges (ex: [e1p1 e1p2 e2p1 e2p2 ...])."))
        , d_surface_triangles(initData(&d_surface_triangles,
                "surface_triangles",
                "List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...])."))
        , d_positions(initData(&d_positions, SofaVecCoord(),
                "position",
                "Position vector of nodes contained in the sparse grid."))
        , d_quads(initData(&d_quads,
                "quads",
                "List of quads contained in the sparse grid (ex: [q1p1 q1p2 q1p3 q1p4 q2p1 ... qnp3 qnp4])."))
        , d_hexahedrons(initData(&d_hexahedrons,
                "hexahedrons",
                "List of hexahedrons contained in the sparse grid (ex: [h1p1 h1p2 h1p3 h1p4 h1p5 ... hnp6 hnp7])."))
{
            addAlias(&d_positions,   "positions");
            addAlias(&d_hexahedrons, "hexahedra");
}

template <typename DataTypes>
void FictitiousGrid<DataTypes>::init() {
    using namespace sofa::core::topology;
    using namespace sofa::core::behavior;

    const auto & triangles = d_surface_triangles.getValue();
    const auto & edges = d_surface_edges.getValue();

    // Let's first check if an iso-surface is specified
    if (d_iso_surface.get()) {
        if (not triangles.empty() or not edges.empty()) {
            msg_warning() << "Both an iso-surface and a surface tesselation (triangles or edges) were specified."
                             "The iso-surface will be used.";
        }
    } else {
        // Let's check if an iso-surface can be found in the current context
        auto iso_surfaces = this->getContext()->template getObjects<IsoSurface<DataTypes>>(BaseContext::Local);
        if (iso_surfaces.size() > 0) {
            if (iso_surfaces.size() > 1) {
                msg_warning() << "Multiple iso-surfaces were found in the context node. The first one, '"
                              << iso_surfaces[0]->getName()
                              << "' will be used. If this is not the one that was needed, please set the parameter '"
                              << d_iso_surface->getPathName()
                              << "' with the path to the correct iso-surface to use.";
            }
            d_iso_surface.set(iso_surfaces[0]);
            msg_info() << "Automatically found the iso-surface '" << d_iso_surface.getPath() << "'.";
        }
    }

    if (not d_iso_surface.get()) {
        BaseMeshTopology * topology = nullptr;
        // Fetch the surface elements (edges in 2D, triangles in 3D)
        if ((Dimension == 2 and edges.empty()) or (Dimension == 3 and triangles.empty())) {
            auto containers = this->getContext()->template getObjects<BaseMeshTopology>(BaseContext::Local);
            auto node = static_cast<const sofa::simulation::Node *> (this->getContext());
            if (containers.empty() || containers.size() > 1) {
                std::vector<BaseMeshTopology *> surface_containers;
                for (auto * c : containers) {
                    if (Dimension == 2 and c->getNbEdges() > 0)
                        surface_containers.emplace_back(c);
                    else if (Dimension == 3 and (c->getNbTriangles() > 0))
                        surface_containers.emplace_back(c);
                }

                if (surface_containers.empty()) {
                    msg_error()
                        << "No surface elements provided and no container were found in the context node '"
                        << node->getPathName() << "'. "
                        << "Please set the parameter '" << d_surface_edges.getName() << "' to the list of triangles"
                        << " if the grid is in 2D, or '" << d_surface_triangles.getName() << "' if the grid is in 3D."
                        << "Note that only triangles are supported for 3D grid.";
                    return;
                } else if (surface_containers.size() > 1) {
                    msg_error()
                        << "Multiple topology containers were found in the context node '" << node->getPathName() << "'."
                        << "Please set the parameter '" << d_surface_edges.getName() << "' to the list of triangles"
                        << " if the grid is in 2D, or '" << d_surface_triangles.getName() << "' if the grid is in 3D."
                        << "Note that only triangles are supported for 3D grid.";
                    return;
                } else {
                    topology = surface_containers[0];
                    if (Dimension == 2) {
                        sofa::helper::WriteAccessor<Data<sofa::type::vector<SofaEdge>>> w_edges = d_surface_edges;
                        for (unsigned int i = 0; i < topology->getNbEdges(); ++i) {
                            w_edges.push_back(topology->getEdge(i));
                        }
                        msg_info() << "Automatically found " << topology->getNbEdges() << " edges in the container '"
                                   << topology->getPathName() << "'.";
                    } else {
                        sofa::helper::WriteAccessor<Data<sofa::type::vector<SofaTriangle >>> w_triangles = d_surface_triangles;
                        for (unsigned int i = 0; i < topology->getNbTriangles(); ++i) {
                            w_triangles.push_back(topology->getTriangle(i));
                        }
                        msg_info() << "Automatically found " << topology->getNbTriangles() << " triangles in the container '"
                                   << topology->getPathName() << "'.";
                    }
                }
            } else {
                topology = containers[0];
                if (Dimension == 2) {
                    if (topology->getNbEdges() > 0) {
                        sofa::helper::WriteAccessor<Data<sofa::type::vector<SofaEdge>>> w_edges = d_surface_edges;
                        for (unsigned int i = 0; i < topology->getNbEdges(); ++i) {
                            w_edges.push_back(topology->getEdge(i));
                        }
                        msg_info() << "Automatically found " << topology->getNbEdges() << " edges in the container '"
                                   << topology->getPathName() << "'.";
                    } else {
                        msg_warning() << "The container '" << topology->getPathName() << "' was found, but contains no edges.";
                    }
                } else {
                    if (topology->getNbTriangles() > 0) {
                        sofa::helper::WriteAccessor<Data<sofa::type::vector<SofaTriangle >>> w_triangles = d_surface_triangles;
                        for (unsigned int i = 0; i < topology->getNbTriangles(); ++i) {
                            w_triangles.push_back(topology->getTriangle(i));
                        }
                        msg_info() << "Automatically found " << topology->getNbTriangles()
                                   << " triangles in the container '"
                                   << topology->getPathName() << "'.";
                    } else {
                        msg_warning() << "The container '" << topology->getPathName() << "' was found, but contains no triangles.";
                    }
                }
            }
        }

        // Fetch the position vector of the surface topology
        if (d_surface_positions.getParent()) {
            // The position vector is a link to another component's vector data field
            sofa::helper::ReadAccessor<Data<SofaVecCoord>> x = d_surface_positions;
            if (x.empty()) {
                msg_error()
                    << "The position vector is linked to '" << d_surface_positions.getLinkPath()
                    << "', but this vector does not contain any node positions.";
                return;
            }
        } else {
            sofa::helper::ReadAccessor<Data<SofaVecCoord>> x = d_surface_positions;
            if (x.empty()) {
                // No position specified, try to find a mechanical object in the context node of the surface topology
                if (d_surface_positions.getParent()) {
                    msg_error() << "The surface boundary positions vector was linked to the field '"
                                << d_surface_positions.getParent()->getLinkPath()
                                << "', but this one does not contain any positions.";
                    return;
                } else {
                    auto mesh_topology = (topology) ? dynamic_cast<sofa::component::topology::container::constant::MeshTopology*>(topology) : nullptr;
                    if (mesh_topology) {
                        d_surface_positions.setParent(&(mesh_topology->seqPoints));
                        msg_info() << "Automatically found the positions vector at '"
                                   << mesh_topology->getPathName() + "." + mesh_topology->seqPoints.getName() << "'.";
                    } else {
                        auto surface_context = (topology) ? topology->getContext() : this->getContext();
                        auto surface_node = dynamic_cast<const sofa::simulation::Node *> (surface_context);
                        auto mechanical_states = surface_context->template getObjects<MechanicalState<DataTypes>>(
                            BaseContext::Local);
                        if (mechanical_states.empty()) {
                            msg_error() << "No positions were set in the '" << d_surface_positions.getName()
                                        << "' data parameter, "
                                        << "and no mechanical objects were found in the context node ('"
                                        << surface_node->getPathName() << "').";
                            return;
                        } else if (mechanical_states.size() > 1) {
                            msg_error() << "No positions were set in the '" << d_surface_positions.getName()
                                        << "' data parameter, "
                                        << "and more than one mechanical objects were found in the context node ('"
                                        << surface_node->getPathName() << "').";
                            return;
                        } else {
                            const MechanicalState<DataTypes> *mo = mechanical_states[0];
                            const auto rest_positions_data = dynamic_cast<Data<SofaVecCoord> *> (mo->findData(
                                "rest_position"));
                            if (!rest_positions_data) {
                                msg_error() << "A mechanical state was found at '" << mo->getPathName()
                                            << "' but does not contain a data field named 'rest_position'.";
                                return;
                            } else {
                                sofa::helper::ReadAccessor<Data<SofaVecCoord>> rest_positions = rest_positions_data;
                                if (rest_positions.empty()) {
                                    msg_error() << "A mechanical state was automatically found at '"
                                                << mo->getPathName()
                                                << "', "
                                                << "but it does not contain any positions in its data field named '"
                                                << rest_positions_data->getName() << "'.";
                                    return;
                                } else {
                                    d_surface_positions.setParent(rest_positions_data);
                                    msg_info() << "Automatically found the positions vector at '"
                                               << mo->getPathName() + "." + rest_positions_data->getName() << "'.";
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // If no bounds are given to the grid, find the bounding box of the surface nodes
    if (not d_min.isSet() or not d_max.isSet()) {
        // Make sure we have a surface mesh and not an implicit field
        if (d_iso_surface.get()) {
            msg_error() << "No bounds (min and max corner) given to the fictitious grid.";
        } else {
            const auto & surface_positions = d_surface_positions.getValue();
            static const SofaFloat max_real = std::numeric_limits<SofaFloat>::max();
            static const SofaFloat min_real = std::numeric_limits<SofaFloat>::lowest();
            SofaVecFloat min, max;
            max.fill(min_real);
            min.fill(max_real);

            for (const auto & p : surface_positions) {
                for (int c=0; c<Dimension; c++) {
                    if (p[c] > max[c])
                        max[c] = static_cast<SofaFloat>(p[c]);
                    if (p[c] < min[c])
                        min[c] = static_cast<SofaFloat>(p[c]);
                }
            }
            d_min.setValue(min);
            d_max.setValue(max);
        }
    }

    // Create the grid
    create_grid();
}


template <typename DataTypes>
void FictitiousGrid<DataTypes>::create_grid()
{
    // Initializing the regular grid
    const auto first_corner = d_min.getValue();
    const auto second_corner = d_max.getValue();
    const auto n = d_n.getValue();

    WorldCoordinates anchor_position;
    Dimensions grid_size;
    for (sofa::Index axis = 0; axis < caribou::geometry::traits<CellElement>::Dimension; ++axis) {
        anchor_position[axis] = std::min(first_corner[axis], second_corner[axis]);
        grid_size[axis] = std::abs(first_corner[axis] - second_corner[axis]);
    }

    Subdivisions grid_n = Subdivisions(&n[0]) - Subdivisions::Ones();

    p_grid = std::make_unique<GridType> (
        anchor_position, grid_n, grid_size
    );

    p_cells_types.resize(p_grid->number_of_cells(), Type::Undefined);
    p_cells.resize(p_grid->number_of_cells());

    // Initialize the full regular grid quadtree (resp. octree) with 0 subdivisions
    for (auto & cell : p_cells) {
        cell.data = std::make_unique<CellData>(Type::Undefined, 1, -1, false);
        cell.childs.reset();
    }

    // 1. Find the cells intersected by the boundary and tag them as "Boundary".
    if (d_iso_surface.get()) {
        tag_intersected_cells_from_implicit_surface();
    } else {
        tag_intersected_cells();
    }

    // 2. Subdivision of the intersected cells
    //    Iteratively subdivide intersected cells until the maximum level of subdivision is reached. Each leaves cells
    //    interesected by the boundary are then tagged as "boundary". Remaining cells are tagged as "Undefined".
    subdivide_intersected_cells();

    // 2. Clustering of cells
    //    Cells of the same type (either undefined or boundary) are grouped in regions where two cells of the same type
    //    but separated by cells of another type cannot be in the same region.
    create_regions_from_same_type_cells();

    // 3. Identifying and tagging of outside regions
    tag_outside_cells();

    // 4. Identifying and tagging remaining regions as inside
    tag_inside_cells();

    // 5. Make sure the grid is valid
    validate_grid();

    // 5. Creating a vector of hexahedral elements containing only inside and boundary cells
    create_sparse_grid();

    // 6. Prepopulating vectors used for display
    populate_drawing_vectors();
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::tag_intersected_cells_from_implicit_surface()
{
    if (!p_grid or p_grid->number_of_nodes() == 0)
        return;

    if (!d_iso_surface.get()) {
        return;
    }

    BEGIN_CLOCK;
    TICK;
    const auto iso_surface = d_iso_surface.get();
    const auto number_of_nodes = p_grid->number_of_nodes();
    const auto number_of_cells = p_grid->number_of_cells();
    std::vector<Type> node_types (number_of_nodes, Type::Undefined);

    // We first compute the type of every nodes
    for (UNSIGNED_INTEGER_TYPE i = 0; i < number_of_nodes; ++i) {
        const auto t = iso_surface->iso_value(p_grid->node(i));
        if (t < 0)
            node_types[i] = Type::Inside;
        else if (t > 0)
            node_types[i] = Type::Outside;
        else
            node_types[i] = Type::Boundary;
    }

    // Once we got the type of the nodes, we compute the type of their cells
    for (UNSIGNED_INTEGER_TYPE cell_index = 0; cell_index < number_of_cells; ++cell_index) {
        const auto node_indices = p_grid->node_indices_of(cell_index);
        UNSIGNED_INTEGER_TYPE types[3] = {0, 0, 0}; // Number of nodes that are Inside, Outside and Boundary, respectively

        for (const auto & node_index : node_indices) {
            types[(UNSIGNED_INTEGER_TYPE) node_types[node_index]]++;
        }

        if (types[(UNSIGNED_INTEGER_TYPE) Type::Inside] == CellElement::NumberOfNodesAtCompileTime) {
            p_cells_types[cell_index] = Type::Inside;
        } else if (types[(UNSIGNED_INTEGER_TYPE) Type::Outside] == CellElement::NumberOfNodesAtCompileTime) {
            p_cells_types[cell_index] = Type::Outside;
        } else {
            p_cells_types[cell_index] = Type::Boundary;
        }
    }

    msg_info() << "Computing the intersections with the surface in "  << std::setprecision(3) << std::fixed
               << TOCK/1000./1000. << " [ms]";
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::create_regions_from_same_type_cells()
{
    BEGIN_CLOCK;
    // At this point, we have all the boundary cells, let's create regions of cells regrouping neighbors cells
    // of the same type. Special attention needs to be spent on boundary cells as if the number of subdivision is
    // greater than zero, we must also add subcells to the regions.
    static constexpr INTEGER_TYPE directions[2] = {-1, 1};

    TICK;

    // First, add all leaf cells to a queue.
    std::queue<Cell*> leaf_cells;
    for (std::size_t i = 0; i < p_grid->number_of_cells(); ++i) {
        for (const auto & lc : get_leaf_cells(&p_cells[i])) {
            leaf_cells.emplace(lc);
        }
    }

    // Iterate over every cells, and for each one, if it isn't yet classified,
    // start a clustering algorithm to fill the region's cells
    p_regions.clear();
    while (not leaf_cells.empty()) {
        Cell * c = leaf_cells.front();
        leaf_cells.pop();

        // If this cell was previously classified, skip it
        if (c->data->region_id > -1) {
            continue;
        }

        // Create a new region
        p_regions.push_back(Region {c->data->type, std::vector<Cell*> ()});
        int current_region_index = static_cast<int>(p_regions.size() - 1);

        // Tag the current cell
        c->data->region_id = current_region_index;
        p_regions[current_region_index].cells.push_back(c);

        // Start the clustering algorithm
        std::queue<Cell *> cluster_cells;
        cluster_cells.push(c);

        while (not cluster_cells.empty()) {
            Cell * cluster_cell = cluster_cells.front();
            cluster_cells.pop();

            // Add the neighbors that aren't already classified
            // For each axis (x, y [, z])
            for (std::size_t axis = 0; axis < Dimension; ++axis) {
                // For each direction (-1 or +1)
                for (const auto &direction : directions) {
                    auto face_neighbors = get_neighbors(cluster_cell, axis, direction);
                    if (face_neighbors.size() == 0) {
                        // The cell is on the grid's boundaries
                        cluster_cell->data->boundary_of_region = true;
                    } else {
                        for (Cell * neighbor_cell : face_neighbors) {
                            if (neighbor_cell->data->type != p_regions[current_region_index].type) {
                                // The cell is on the boundary if its region
                                cluster_cell->data->boundary_of_region = true;
                            } else {
                                // They are of the same type, hence they are in the same region
                                // Let's make sure of it (maybe this should be on a debug assertion instead)
                                if (neighbor_cell->data->region_id != -1 and neighbor_cell->data->region_id != static_cast<int>(current_region_index)) {
                                    throw std::runtime_error("Two neighboring leaf-cells are of the same type, but not of the same region. This should not happen.");
                                } else if (neighbor_cell->data->region_id == -1) {
                                    // Tag the neighbor cell and add it to the list so we can tag its neighbors.
                                    neighbor_cell->data->region_id = current_region_index;
                                    p_regions[current_region_index].cells.emplace_back(neighbor_cell);
                                    cluster_cells.emplace(neighbor_cell);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    msg_info() << "Computing the " << p_regions.size() << " cells regions in " << std::fixed << std::setprecision(3)
               << TOCK / 1000. / 1000.
               << " [ms]";
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::tag_outside_cells()
{
    BEGIN_CLOCK;
    // At this point, all cells have been classified into distinct regions of either boundary type or undefined type.
    // The regions of undefined type which are surrounded by the grid's boundaries are tagged as outside cells

    TICK;
    const auto & upper_grid_boundary = p_grid->N();
    for (auto & region : p_regions) {
        if (region.type != Type::Undefined) {
            continue;
        }

        for (Cell *cell : region.cells) {
            auto is_boundary_of_face = std::bitset<6>().flip();

            // First, check if one of the subcell's face is the boundary of all its parent cells
            CellIndex index = cell->index;
            Cell *p = cell;
            while (p->parent) {
                const auto &coordinates = subcell_coordinates[index];
                for (UNSIGNED_INTEGER_TYPE axis = 0; axis < Dimension; ++axis) {
                    is_boundary_of_face[2 * axis + 0] = is_boundary_of_face[2 * axis + 0] &&  (coordinates[axis] == 0);
                    is_boundary_of_face[2 * axis + 1] = is_boundary_of_face[2 * axis + 1] &&  (coordinates[axis] == 1);
                }
                p = p->parent;
                index = p->index;
            }

            // Next, check if any of the top parent's faces are at the boundary of the grid
            const auto coordinates = p_grid->cell_coordinates_at(index);
            for (UNSIGNED_INTEGER_TYPE axis = 0; axis < Dimension; ++axis) {
                is_boundary_of_face[2 * axis + 0] = is_boundary_of_face[2 * axis + 0] && (coordinates[axis] == 0);
                is_boundary_of_face[2 * axis + 1] = is_boundary_of_face[2 * axis + 1] &&
                                                    ((unsigned) coordinates[axis] == upper_grid_boundary[axis]-1);
            }

            // If the subcell has a face that is both on the boundary of all of its parents, and is also boundary
            // of the grid, than this subcell is outside of the surface and we can tag all cells in the group as outside
            // cells
            if (is_boundary_of_face.any()) {
                region.type = Type::Outside;
                for (Cell *c : region.cells) {
                    c->data->type = Type::Outside;
                }
                break;
            }
        }
    }
    msg_info() << "Computing the outside regions types in " << std::setprecision(3)
               << TOCK / 1000. / 1000.
               << " [ms]";
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::tag_inside_cells()
{
    // Finally, we have both outside and boundary regions, the remaining regions are inside.
    BEGIN_CLOCK;
    TICK;
    for (auto & region : p_regions) {
        if (region.type != Type::Undefined) {
            continue;
        }

        region.type = Type::Inside;
        for (const auto & cell : region.cells) {
            cell->data->type = Type::Inside;
        }
    }
    msg_info() << "Computing the inside regions types in " << std::setprecision(3)
               << TOCK / 1000. / 1000.
               << " [ms]";
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::validate_grid() {
    UNSIGNED_INTEGER_TYPE number_of_inside_regions = 0,
                          number_of_outside_regions = 0,
                          number_of_boundary_regions = 0,
                          number_of_undefined_regions = 0;

    FLOATING_POINT_TYPE volume = 0.;
    FLOATING_POINT_TYPE cell_volume = p_grid->H().prod();
    for (const Region & region : p_regions) {
        if (region.type == Type::Inside) {
            ++number_of_inside_regions;
        } else if (region.type == Type::Outside) {
            ++number_of_outside_regions;
        } else if (region.type == Type::Boundary) {
            ++number_of_boundary_regions;
        } else {
            ++number_of_undefined_regions;
        }
        for (const Cell * cell : region.cells) {
            if (not cell->is_leaf()) {
                throw std::runtime_error("Some cells of a region aren't leaf cells. This should not happen.");
            }
            volume += cell->data->weight * cell_volume;
        }
    }

    FLOATING_POINT_TYPE grid_volume = p_grid->size().prod();
    if (std::abs(volume - grid_volume) > grid_volume*1e-8) {
        std::stringstream ss;
        ss << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
           << "The grid volume (" << grid_volume << ") does not match the volume of the region leaf-cells ("
           << volume << "). This should not happen.";
        throw std::runtime_error(ss.str());
    }

    msg_info() << "There are " << number_of_inside_regions << " inside regions";
    msg_info() << "There are " << number_of_outside_regions << " outside regions";
    msg_info() << "There are " << number_of_boundary_regions << " boundary regions";
    msg_info() << "There are " << number_of_undefined_regions << " undefined regions";

    if (number_of_inside_regions == 0) {
        msg_warning() << "There are no cells lying completely inside the surface. You might want to refine the mesh "
                      << "or make sure that there are no holes in the surface mesh.";
    }


}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::create_sparse_grid()
{
    BEGIN_CLOCK;
    TICK;

    const auto & volume_threshold = d_volume_threshold.getValue();
    std::vector<bool> use_cell(p_grid->number_of_cells(), false);
    std::vector<bool> use_node(p_grid->number_of_nodes(), false);

    std::map<UNSIGNED_INTEGER_TYPE, UNSIGNED_INTEGER_TYPE> volume_ratios;
    UNSIGNED_INTEGER_TYPE ignored_cells_count = 0;
    FLOATING_POINT_TYPE real_volume = 0.;
    FLOATING_POINT_TYPE cell_volume = CellElement::NumberOfGaussNodesAtCompileTime*p_grid->cell_at(0).jacobian(CellElement::LocalCoordinates::Zero()).determinant();

    // 1. Locate all cells that are within the surface boundaries and their nodes.
    for (UNSIGNED_INTEGER_TYPE cell_id = 0; cell_id < p_grid->number_of_cells(); ++cell_id) {
        const Cell & cell = p_cells[cell_id];
        if (not cell.is_leaf() or cell.data->type != Type::Outside) {

            const FLOATING_POINT_TYPE weight = get_cell_weight(cell);
            real_volume += cell_volume*weight;

            const auto ratio = static_cast<UNSIGNED_INTEGER_TYPE> (weight*10)*10;
            if (volume_ratios.find(ratio) == volume_ratios.end())
                volume_ratios[ratio] = 1;
            else
                volume_ratios[ratio] += 1;

            if (weight >= volume_threshold) {
                use_cell[cell_id] = true;
                for (const auto node_index : p_grid->node_indices_of(cell_id)) {
                    use_node[node_index] = true;
                }
            } else {
                ++ignored_cells_count;
            }
        }
    }


    // 2. Add the sparse nodes and create the bijection between sparse nodes and full grid nodes
    sofa::helper::WriteAccessor<Data< SofaVecCoord >> positions = d_positions;
    positions.clear();
    positions.reserve(p_grid->number_of_nodes());

    p_node_index_in_sparse_grid.clear();
    p_node_index_in_sparse_grid.resize(p_grid->number_of_nodes(), -1);
    p_node_index_in_grid.clear();
    p_node_index_in_grid.reserve(p_grid->number_of_nodes());

    for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < p_grid->number_of_nodes(); ++node_id) {
        if (use_node[node_id]) {
            const auto position = p_grid->node(node_id);
            if constexpr (Dimension == 2) {
                positions.wref().emplace_back(position[0], position[1]);
            } else {
                positions.wref().emplace_back(position[0], position[1], position[2]);
            }
            p_node_index_in_grid.emplace_back(node_id);
            p_node_index_in_sparse_grid[node_id] = (signed) positions.size() - 1;
        }
    }

    // 3. Add the sparse cells and create the bijection between sparse cells and full grid cells
    sofa::helper::WriteAccessor<Data < sofa::type::vector<SofaHexahedron> >> hexahedrons = d_hexahedrons;
    sofa::helper::WriteAccessor<Data < sofa::type::vector<SofaQuad > >> quads = d_quads;
    if (Dimension == 2) {
        quads.clear();
        quads.wref().reserve(p_grid->number_of_cells());
    } else {
        hexahedrons.clear();
        hexahedrons.wref().reserve(p_grid->number_of_cells());
    }

    p_cell_index_in_sparse_grid.clear();
    p_cell_index_in_sparse_grid.resize(p_grid->number_of_cells(), -1);
    p_cell_index_in_grid.clear();
    p_cell_index_in_grid.reserve(p_grid->number_of_cells());

    for (UNSIGNED_INTEGER_TYPE cell_id = 0; cell_id < p_grid->number_of_cells(); ++cell_id) {
        if (not use_cell[cell_id]) {
            continue;
        }

        const auto node_indices = p_grid->node_indices_of(cell_id);
        if (Dimension == 2) {
            quads.wref().emplace_back(
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[0]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[1]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[2]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[3]])
            );
            p_cell_index_in_sparse_grid[cell_id] = (signed) quads.size() - 1;
        } else {
            hexahedrons.wref().emplace_back(
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[0]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[1]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[2]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[3]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[4]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[5]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[6]]),
                static_cast<sofa::Index>(p_node_index_in_sparse_grid[node_indices[7]])
            );
            p_cell_index_in_sparse_grid[cell_id] = (signed) hexahedrons.size() - 1;
        }

        p_cell_index_in_grid.emplace_back(cell_id);
    }

    positions.wref().shrink_to_fit();
    hexahedrons.wref().shrink_to_fit();
    quads.wref().shrink_to_fit();
    p_node_index_in_grid.shrink_to_fit();
    p_cell_index_in_grid.shrink_to_fit();

    msg_info() << "Creating the sparse grid in " << std::setprecision(3)
               << TOCK / 1000. / 1000.
               << " [ms]";

    msg_info() << "Volume of the sparse grid is " << real_volume;
    for (const auto &p : volume_ratios) {
        if (p.first == 100) {
            msg_info() << p.second << " cells with a volume ratio of 100%";
        } else {
            msg_info() << p.second << " cells with a volume ratio between " << p.first << "% and " << p.first + 10 << "%";
        }
    }
    if (ignored_cells_count > 0) {
        msg_info() << ignored_cells_count << " cells with a volume ratio less than " << static_cast<unsigned int>(volume_threshold*100) << "% (IGNORED CELLS)";
    }
}

template <typename DataTypes>
auto
FictitiousGrid<DataTypes>::get_subcells_elements(const CellElement & e) const -> std::array<CellElement, (unsigned) 1 << Dimension>
{
    const Dimensions h = e.size() /  2;
    if constexpr (Dimension == 2) {
        return {{
                    CellElement(e.world_coordinates(LocalCoordinates {-0.5, -0.5}), h),
                    CellElement(e.world_coordinates(LocalCoordinates {+0.5, -0.5}), h),
                    CellElement(e.world_coordinates(LocalCoordinates {-0.5, +0.5}), h),
                    CellElement(e.world_coordinates(LocalCoordinates {+0.5, +0.5}), h)
        }};
    } else {
        return {{
                   CellElement(e.world_coordinates(LocalCoordinates {-0.5, -0.5, -0.5}), h),
                   CellElement(e.world_coordinates(LocalCoordinates {+0.5, -0.5, -0.5}), h),
                   CellElement(e.world_coordinates(LocalCoordinates {-0.5, +0.5, -0.5}), h),
                   CellElement(e.world_coordinates(LocalCoordinates {+0.5, +0.5, -0.5}), h),
                   CellElement(e.world_coordinates(LocalCoordinates {-0.5, -0.5, +0.5}), h),
                   CellElement(e.world_coordinates(LocalCoordinates {+0.5, -0.5, +0.5}), h),
                   CellElement(e.world_coordinates(LocalCoordinates {-0.5, +0.5, +0.5}), h),
                   CellElement(e.world_coordinates(LocalCoordinates {+0.5, +0.5, +0.5}), h)
        }};
    }
}

template <typename DataTypes>
std::vector<typename FictitiousGrid<DataTypes>::Cell *>
FictitiousGrid<DataTypes>::get_leaf_cells(const Cell * cell) const
{
    const auto & number_of_subdivision = d_number_of_subdivision.getValue();
    std::vector<typename FictitiousGrid<DataTypes>::Cell *> leafs;

    // Reserve the space for a maximum of 2^(nbSudiv*2) leafs in 2D, 2^(nbSudiv*3)  leafs in 3D
    leafs.reserve(static_cast<UNSIGNED_INTEGER_TYPE>(1) << (number_of_subdivision*Dimension));

    std::queue<const Cell *> cells;
    cells.emplace(cell);

    while (not cells.empty()) {
        const Cell * c = cells.front();
        cells.pop();
        if (not c->childs) {
            leafs.emplace_back(const_cast<Cell *>(c));
        } else {
            for (const Cell & child : *(c->childs)) {
                cells.emplace(&child);
            }
        }
    }

    leafs.shrink_to_fit();

    return leafs;
}

template <typename DataTypes>
std::vector<typename FictitiousGrid<DataTypes>::Cell *>
FictitiousGrid<DataTypes>::get_neighbors(const Cell * cell) const
{
    std::vector<Cell *> neighbors;

    static constexpr INTEGER_TYPE directions[2] = {-1, 1};


    // For each axis (x, y [, z])
    for (std::size_t axis = 0; axis < Dimension; ++axis) {
        // For each direction (-1 or +1)
        for (const auto & direction : directions) {
            const auto face_neighbors = get_neighbors(cell, axis, direction);
            neighbors.insert(neighbors.end(), face_neighbors.begin(), face_neighbors.end());
        }
    }

    return neighbors;
}

template <typename DataTypes>
std::vector<typename FictitiousGrid<DataTypes>::Cell *>
FictitiousGrid<DataTypes>::get_neighbors(const Cell * cell, UNSIGNED_INTEGER_TYPE axis, INTEGER_TYPE direction) const
{
    const auto & upper_grid_boundary = p_grid->N();
    std::vector<Cell *> neighbors;

    const Cell * p = cell;
    CellIndex index;
    bool found_suitable_parent = false;
    std::stack<GridCoordinates> path;


    // Go up in the quadtree (resp. octree) until we can move in the axis-direction, and save the path
    // while doing it so we can follow it back when going down in the neighbor cell
    do {
        index = p->index;
        if (p->parent) {
            const GridCoordinates &coordinates = subcell_coordinates[index];
            if (IN_CLOSED_INTERVAL(0, coordinates[axis] + direction, 1)) {
                found_suitable_parent = true;
            }
            path.emplace(coordinates);
        }
        p = p->parent;
    } while (p and not found_suitable_parent);

    // If we can no longer move up, we are at the cell top most parent. Check if we can move in the
    // axis-direction within the grid
    if (not found_suitable_parent) {
        const GridCoordinates coordinates = p_grid->cell_coordinates_at(index);
        const int new_coordinate = static_cast<int>(coordinates[axis] + direction);
        const int upper_limit    = static_cast<int>(upper_grid_boundary[axis] - 1);
        if (0 <= new_coordinate and new_coordinate <= upper_limit) {
            // We got a neighbor cell within the grid's limit
            auto new_coordinates = coordinates;
            new_coordinates[axis] += direction;
            const Cell & c = p_cells[p_grid->cell_index_at(new_coordinates)];
            p = &c;
        }
    }

    if (p) {
        // We found a neighbors cell within the limits of the cell quadtree (resp. octree) or the grid.
        // Let's go down following the reverse path until we get to the same level as the queried cell
        while (not path.empty() and p->childs) {
            GridCoordinates coordinates = path.top();
            // We reverse the axis coordinate so that we can get cells at the opposite of the queried one
            const auto new_coordinate = (coordinates[axis] - direction) % 2;
            coordinates[axis] = new_coordinate*new_coordinate;
            index = coordinates[0] + coordinates[1]*2;
            if constexpr (Dimension == 3) {
                index += coordinates[2]*2*2;
            }
            p = &(*p->childs)[index];
            path.pop();
        }

        // We are now at the same level of the queried cells, on the opposite side in the axis direction.
        // If this cell still has children, lets find those that are neighbors to the queried cell.
        if (p->childs) {

            // Let's set the axis direction so that we only get children facing the queried cell.
            INTEGER_TYPE search_axis = (direction < 0) ? 1 : 0;

            std::stack<Cell*> children;
            for (auto & child : *p->childs) {
                children.emplace(&child);
            }
            while (not children.empty()) {
                const Cell * c = children.top();
                children.pop();
                GridCoordinates coordinates = subcell_coordinates[c->index];
                if (coordinates[axis] == search_axis) {
                    if (c->data) {
                        neighbors.emplace_back(const_cast<Cell *>(c));
                    } else {
                        for (auto & child : *c->childs) {
                            children.emplace(&child);
                        }
                    }
                }
            }
        } else {
            // No further childs, let's add the cell as a neighbor
            neighbors.emplace_back(const_cast<Cell *>(p));
        }
    }
    return neighbors;
}

template <typename DataTypes>
inline FLOATING_POINT_TYPE FictitiousGrid<DataTypes>::get_cell_weight(const Cell & cell) const
{
    if (cell.is_leaf()) {
        if (cell.data->type == Type::Inside or cell.data->type == Type::Boundary)
            return cell.data->weight;
        else
            return 0.;
    } else {
        FLOATING_POINT_TYPE w = 0;
        for (const auto & c : *(cell.childs)) {
            w += get_cell_weight(c);
        }
        return w;
    }
}

template <typename DataTypes>
std::vector<std::pair<typename FictitiousGrid<DataTypes>::LocalCoordinates, FLOATING_POINT_TYPE>>
FictitiousGrid<DataTypes>::get_gauss_nodes_of_cell(const CellIndex & sparse_cell_index) const
{
    return get_gauss_nodes_of_cell(sparse_cell_index, d_number_of_subdivision.getValue());
}

template <typename DataTypes>
std::vector<std::pair<typename FictitiousGrid<DataTypes>::LocalCoordinates, FLOATING_POINT_TYPE>>
FictitiousGrid<DataTypes>::get_gauss_nodes_of_cell(const CellIndex & sparse_cell_index, const UNSIGNED_INTEGER_TYPE maximum_level) const
{
    using Level = UNSIGNED_INTEGER_TYPE;

    const auto cell_index = p_cell_index_in_grid[sparse_cell_index];

    const auto get_gauss_cells = [this, maximum_level] (const CellElement & e, const Cell & c, const Level l) {
        auto get_gauss_cells_impl = [this, maximum_level] (const CellElement & e, const Cell & c, const Level l, auto & ref) -> std::vector<std::pair<LocalCoordinates, FLOATING_POINT_TYPE>> {
            std::vector<std::pair<LocalCoordinates, FLOATING_POINT_TYPE>> gauss_nodes;
            // Reserve the space for a maximum of 2^(level*2) nodes in 2D, 2^(level*3) nodes in 3D
            gauss_nodes.reserve((UNSIGNED_INTEGER_TYPE) 1 << (((maximum_level-l)+1)*Dimension));

            if (c.is_leaf()) {
                for (const auto & gauss_node : e.gauss_nodes() ) {
                    if (c.data->type == Type::Inside or c.data->type == Type::Boundary)
                        gauss_nodes.emplace_back(gauss_node.position, gauss_node.weight*c.data->weight);
                    else
                        gauss_nodes.emplace_back(gauss_node.position, 0);
                }
            } else {
                if (l >= maximum_level) {
                    const auto & child_cells = *(c.childs);
                    for (UNSIGNED_INTEGER_TYPE i = 0; i < CellElement::NumberOfGaussNodesAtCompileTime; ++i) {
                        const auto gauss_node = e.gauss_node(i);
                        const auto weight = CellElement::NumberOfGaussNodesAtCompileTime*get_cell_weight(child_cells[i]);
                        gauss_nodes.emplace_back(LocalCoordinates(gauss_node.position), gauss_node.weight*weight);
                    }
                } else {
                    const auto child_elements = get_subcells_elements(e);
                    const auto & child_cells = *(c.childs);
                    for (UNSIGNED_INTEGER_TYPE i = 0; i < child_cells.size(); ++i) {
                        const auto child_gauss_nodes = ref(child_elements[i], child_cells[i], l+1, ref);
                        for (const auto & child_gauss_node : child_gauss_nodes) {
                            gauss_nodes.emplace_back(child_elements[i].world_coordinates(child_gauss_node.first), child_gauss_node.second);
                        }
                    }
                }
            }
            gauss_nodes.shrink_to_fit();
            return gauss_nodes;
        };
        return get_gauss_cells_impl(e, c, l, get_gauss_cells_impl);
    };

    const CellElement top_element = p_grid->cell_at(cell_index);
    const FLOATING_POINT_TYPE detJ = top_element.jacobian(CellElement::LocalCoordinates::Zero()).determinant();
    std::vector<std::pair<LocalCoordinates, FLOATING_POINT_TYPE>> gauss_nodes = get_gauss_cells(CellElement(), p_cells[cell_index], 0);

    for (auto & n : gauss_nodes) {
        n.second *= detJ;
    }

    return gauss_nodes;
}

template <typename DataTypes>
auto FictitiousGrid<DataTypes>::get_type_at(const WorldCoordinates & p) const -> Type {
    if (p_grid->contains(p)) {
        const auto cells = p_grid->cells_around(p);

        if (cells.size() == 1) {
            return p_cells_types[cells[0]];
        }

        // The position p is on the boundary of multiple cells, gather the different cells types
        INTEGER_TYPE types = 0;
        for (const auto & cell_index : cells) {
            types |= static_cast<INTEGER_TYPE>(p_cells_types[cell_index]);
        }

        if (types & static_cast<INTEGER_TYPE>(Type::Boundary)) {
            // If one of the cells around p is of type boundary, return the type boundary
            return Type::Boundary;
        } else if((types & static_cast<INTEGER_TYPE>(Type::Inside)) and (types & static_cast<INTEGER_TYPE>(Type::Outside))) {
            // If one of the cells around p is of type inside, and another one is of type outside, return the type boundary
            return Type::Boundary;
        } else {
            // We cannot decide which type it is...this should never happen. Report it.
            std::ostringstream error;
            if (cells.empty()) {
                // Normally should have been caught by the grid->contains test
                error << "The position " << p << " was found within the boundaries of the grid, but is not part of any grid cells.";
            } else {
                error << "The position " << p << " is contained inside multiple cells of different types, and"
                                                 " the type of the position cannot be determined.";
            }

            throw std::runtime_error(error.str());
        }
    } else {
        return Type::Outside;
    }
}

template <typename DataTypes>
auto FictitiousGrid<DataTypes>::cell_volume_ratio_distribution(UNSIGNED_INTEGER_TYPE number_of_decimals) const -> std::map<FLOATING_POINT_TYPE, std::vector<CellIndex>> {
    std::map<FLOATING_POINT_TYPE, std::vector<CellIndex>> volume_ratios;
    const auto number_of_cells = p_grid->number_of_cells();

    const FLOATING_POINT_TYPE step = std::pow(10, number_of_decimals);

    for (CellIndex cell_id = 0; cell_id < static_cast<CellIndex>(number_of_cells); ++cell_id) {
        const auto & cell = p_cells[cell_id];
        FLOATING_POINT_TYPE ratio;

        switch (p_cells_types[cell_id]) {
            case Type::Outside:
                ratio = 0.;
                break;
            case Type::Inside:
                ratio = 1.;
                break;
            case Type::Boundary:
                ratio = (number_of_decimals == 0)
                        ? get_cell_weight(cell)
                        : std::round(get_cell_weight(cell)*step) / step;
                break;
            case Type::Undefined:
            default:
                continue;
        }

        volume_ratios[ratio].emplace_back(cell_id);
    }
    return volume_ratios;
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::populate_drawing_vectors()
{
    BEGIN_CLOCK;
    p_drawing_nodes_vector.resize(p_grid->number_of_nodes());
    p_drawing_edges_vector.resize(p_grid->number_of_edges()*2);
    // Reset the vector in case we got multiple initializations
    p_drawing_subdivided_edges_vector.clear();
    p_drawing_cells_vector.clear();
    p_drawing_subdivided_edges_vector.resize(p_regions.size());
    p_drawing_cells_vector.resize(p_regions.size());
    const double & scale = d_draw_scale.getValue();

    TICK;
    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_nodes(); ++i) {
        const auto p = p_grid->node(i);
        if (Dimension == 2) {
            p_drawing_nodes_vector[i] = {p[0], p[1], 0};
        } else {
            p_drawing_nodes_vector[i] = {p[0], p[1], p[2]};
        }
    }

    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_edges(); ++i) {
        const auto edge = p_grid->edge(i);
        const auto p0 = p_grid->node(edge[0]);
        const auto p1 = p_grid->node(edge[1]);
        if (Dimension == 2) {
            p_drawing_edges_vector[i*2] = {p0[0], p0[1], 0};
            p_drawing_edges_vector[i*2+1] = {p1[0], p1[1], 0};
        } else {
            p_drawing_edges_vector[i*2] = {p0[0], p0[1], p0[2]};
            p_drawing_edges_vector[i*2+1] = {p1[0], p1[1], p1[2]};
        }
    }

    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_cells(); ++i) {
        std::queue<std::pair<const Cell *, CellElement >> cells;
        cells.emplace(&p_cells[i], p_grid->cell_at(i));

        while (not cells.empty()) {
            const Cell * c = std::get<0>(cells.front());
            const CellElement & e = std::get<1>(cells.front());

            if (not c->is_leaf()) {
                const auto & childs_elements = get_subcells_elements(e);
                for (UNSIGNED_INTEGER_TYPE j = 0; j < (*c->childs).size(); ++j) {
                    cells.emplace(&(*c->childs)[j], childs_elements[j]);
                }
            } else if (c->data->boundary_of_region) {
                const auto & region_id = c->data->region_id;
                const WorldCoordinates center = e.center();
                for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < ((unsigned) 1 << Dimension); ++node_id) {
                    const WorldCoordinates p = (center + (e.node(node_id) - center)*scale).transpose();
                    if (Dimension == 2) {
                        p_drawing_cells_vector[region_id].emplace_back(p[0], p[1], 0);
                    } else {
                        p_drawing_cells_vector[region_id].emplace_back(p[0], p[1], p[2]);
                    }
                }
                for (const auto &edge : e.edges()) {
                    const WorldCoordinates p0 = (center + (e.node(edge[0]) - center)*scale).transpose();
                    const WorldCoordinates p1 = (center + (e.node(edge[1]) - center)*scale).transpose();
                    if (Dimension == 2) {
                        p_drawing_subdivided_edges_vector[region_id].emplace_back(p0[0], p0[1], 0);
                        p_drawing_subdivided_edges_vector[region_id].emplace_back(p1[0], p1[1], 0);
                    } else {
                        p_drawing_subdivided_edges_vector[region_id].emplace_back(p0[0], p0[1], p0[2]);
                        p_drawing_subdivided_edges_vector[region_id].emplace_back(p1[0], p1[1], p1[2]);
                    }
                }
            }

            cells.pop();
        }
    }
    msg_info() << "Populating the drawing vectors in " << std::setprecision(3) << std::fixed
               << TOCK / 1000. / 1000.
               << " [ms]";
}

static const unsigned long long kelly_colors_hex[] = {
    0xFFFFB300, // Vivid Yellow
    0xFF803E75, // Strong Purple
    0xFFFF6800, // Vivid Orange
    0xFFA6BDD7, // Very Light Blue
    0xFFC10020, // Vivid Red
    0xFFCEA262, // Grayish Yellow
    0xFF817066, // Medium Gray

    // The following don't work well for people with defective color vision
    0xFF007D34, // Vivid Green
    0xFFF6768E, // Strong Purplish Pink
    0xFF00538A, // Strong Blue
    0xFFFF7A5C, // Strong Yellowish Pink
    0xFF53377A, // Strong Violet
    0xFFFF8E00, // Vivid Orange Yellow
    0xFFB32851, // Strong Purplish Red
    0xFFF4C800, // Vivid Greenish Yellow
    0xFF7F180D, // Strong Reddish Brown
    0xFF93AA00, // Vivid Yellowish Green
    0xFF593315, // Deep Yellowish Brown
    0xFFF13A13, // Vivid Reddish Orange
    0xFF232C16, // Dark Olive Green
};

template <typename DataTypes>
void FictitiousGrid<DataTypes>::draw(const sofa::core::visual::VisualParams* vparams)
{
    using Color = sofa::type::RGBAColor;
    if (!p_grid)
        return;

    const auto & draw_boundary_cells = d_draw_boundary_cells.getValue();
    const auto & draw_outside_cells = d_draw_outside_cells.getValue();
    const auto & draw_inside_cells = d_draw_inside_cells.getValue();

    vparams->drawTool()->disableLighting();


//    vparams->drawTool()->drawPoints(p_drawing_nodes_vector, 1, Color(1, 0, 0, 1));
//    vparams->drawTool()->drawLines(p_drawing_edges_vector, 2, Color(1,0,0, 1));
//    vparams->drawTool()->drawLines(p_drawing_subdivided_edges_vector, 0.5, Color(1,1,1, 1));

    for (UNSIGNED_INTEGER_TYPE region_id = 0; region_id < p_drawing_cells_vector.size(); ++region_id) {
        if ((p_regions[region_id].type == Type::Boundary and draw_boundary_cells) or
            (p_regions[region_id].type == Type::Outside and draw_outside_cells) or
            (p_regions[region_id].type == Type::Inside and draw_inside_cells)) {

            const auto &hex_color = kelly_colors_hex[region_id % 20];
            const Color cell_color(
                (float) (((unsigned char) (hex_color >> (unsigned) 16)) / 255.),
                (float) (((unsigned char) (hex_color >> (unsigned) 8)) / 255.),
                (float) (((unsigned char) (hex_color >> (unsigned) 0)) / 255.),
                (float) (1));
            const Color edge_color(0, 0, 0, 1);

            if (! vparams->displayFlags().getShowWireFrame()) {
                if (Dimension == 2) {
                    vparams->drawTool()->drawQuads(p_drawing_cells_vector[region_id], cell_color);
                } else {
                    vparams->drawTool()->drawHexahedra(p_drawing_cells_vector[region_id], cell_color);
                }
            }
            vparams->drawTool()->drawLines(p_drawing_subdivided_edges_vector[region_id], 0.5, edge_color);
        }
    }

    vparams->drawTool()->restoreLastState();
}

} // namespace SofaCaribou::topology
