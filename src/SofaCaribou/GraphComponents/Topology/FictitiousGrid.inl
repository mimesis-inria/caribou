#ifndef SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL
#define SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL

#include <SofaCaribou/GraphComponents/Topology/FictitiousGrid.h>

#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseTopology/MeshTopology.h>
#include <sofa/simulation/Node.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <stack>
#include <queue>
#include <iomanip>
#include <chrono>

#define BEGIN_CLOCK ;std::chrono::steady_clock::time_point __time_point_begin;
#define TICK ;__time_point_begin = std::chrono::steady_clock::now();
#define TOCK (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - __time_point_begin).count())


namespace SofaCaribou::GraphComponents::topology {

template <typename DataTypes>
FictitiousGrid<DataTypes>::FictitiousGrid()
        : d_n(initData(&d_n, SofaVecInt(),"n","Grid resolution [nx, ny, nz] where nx, ny or nz are the number of nodes in the x,y and z directions."))
        , d_min(initData(&d_min, SofaVecFloat(),"min","First corner node position of the grid's bounding box."))
        , d_max(initData(&d_max, SofaVecFloat(),"max","Second corner node position of the grid's bounding box."))
        , d_number_of_subdivision(initData(&d_number_of_subdivision,
                "maximum_number_of_subdivision_levels",
                0,
                "Number of subdivision levels of the boundary cells (one level split the cell in 4 subcells in 2D, and 8 subcells in 3D)."))
        , d_use_implicit_surface(initData(&d_use_implicit_surface,
                bool(false),
                "use_implicit_surface",
                "Use an implicit surface instead of a tessellated surface. If true, the callback function is_inside must be defined."))
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
                "positions",
                "Position vector of nodes contained in the sparse grid."))
        , d_hexahedrons(initData(&d_hexahedrons,
                "hexahedrons",
                "List of hexahedrons contained in the sparse grid (ex: [h1p1 h1p2 h1p3 h1p4 h1p5 ... hnp6 hnp7])."))
{
}

template <typename DataTypes>
void FictitiousGrid<DataTypes>::init() {
    using namespace sofa::core::topology;
    using namespace sofa::core::behavior;

    const auto & triangles = d_surface_triangles.getValue();
    const auto & edges = d_surface_edges.getValue();

    if (not d_use_implicit_surface.getValue()) {
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
                        sofa::helper::WriteAccessor<Data<sofa::helper::vector<SofaEdge>>> w_edges = d_surface_edges;
                        for (unsigned int i = 0; i < topology->getNbEdges(); ++i) {
                            w_edges.push_back(topology->getEdge(i));
                        }
                        msg_info() << "Automatically found " << topology->getNbEdges() << " edges in the container '"
                                   << topology->getPathName() << "'.";
                    } else {
                        sofa::helper::WriteAccessor<Data<sofa::helper::vector<SofaTriangle >>> w_triangles = d_surface_triangles;
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
                        sofa::helper::WriteAccessor<Data<sofa::helper::vector<SofaEdge>>> w_edges = d_surface_edges;
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
                        sofa::helper::WriteAccessor<Data<sofa::helper::vector<SofaTriangle >>> w_triangles = d_surface_triangles;
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
                    auto mesh_topology = (topology) ? dynamic_cast<sofa::component::topology::MeshTopology*>(topology) : nullptr;
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

        if (not d_min.isSet() or not d_max.isSet()) {
            sofa::helper::ReadAccessor<Data<SofaVecCoord>> surface_positions = d_surface_positions;
            const auto bbox = compute_bbox_from(surface_positions.ref());

            d_min.setValue(bbox.first);
            d_max.setValue(bbox.second);

            msg_info() << "The size of the grid was automatically computed from the bounding box of the surface positions. "
                       << "If you wish to change this value, use the '" << d_min.getName() << "' and '" << d_max.getName() << "' "
                       << "attributes of this component.";
        }
    }

    // Create the grid
    create_grid();
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::compute_cell_types_from_implicit_surface()
{
    if (!p_grid or p_grid->number_of_nodes() == 0)
        return;

    if (!p_implicit_test_callback) {
        return;
    }

    const auto number_of_nodes = p_grid->number_of_nodes();
    const auto number_of_cells = p_grid->number_of_cells();

    // We first compute the type of every nodes
    for (UNSIGNED_INTEGER_TYPE i = 0; i < number_of_nodes; ++i) {
        const auto t = p_implicit_test_callback(p_grid->node(i));
        if (t < 0)
            p_node_types[i] = Type::Inside;
        else if (t > 0)
            p_node_types[i] = Type::Outside;
        else
            p_node_types[i] = Type::Boundary;
    }

    // Once we got the type of the nodes, we compute the type of their cells
    for (UNSIGNED_INTEGER_TYPE cell_index = 0; cell_index < number_of_cells; ++cell_index) {
        const auto node_indices = p_grid->node_indices_of(cell_index);
        UNSIGNED_INTEGER_TYPE types[3] = {0, 0, 0};

        for (const auto & node_index : node_indices) {
            types[(UNSIGNED_INTEGER_TYPE) p_node_types[node_index]]++;
        }

        if (types[(UNSIGNED_INTEGER_TYPE) Type::Inside] == caribou::traits<CellElement>::NumberOfNodes) {
            p_cells_types[cell_index] = Type::Inside;
        } else if (types[(UNSIGNED_INTEGER_TYPE) Type::Outside] == caribou::traits<CellElement>::NumberOfNodes) {
            p_cells_types[cell_index] = Type::Outside;
        } else {
            p_cells_types[cell_index] = Type::Boundary;
        }
    }
}

template <typename DataTypes>
std::vector<typename FictitiousGrid<DataTypes>::Cell *>
FictitiousGrid<DataTypes>::get_neighbors(Cell * cell)
{
    const auto & upper_grid_boundary = p_grid->N();
    std::vector<Cell *> neighbors;

    static constexpr INTEGER_TYPE directions[2] = {-1, 1};


    // For each axis (x, y [, z])
    for (std::size_t axis = 0; axis < Dimension; ++axis) {
        // For each direction (-1 or +1)
        for (const auto & direction : directions) {

            Cell * p = cell;
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
                const GridCoordinates coordinates = p_grid->grid_coordinates_at(index);
                const int new_coordinate = coordinates[axis] + direction;
                const int upper_limit = upper_grid_boundary[axis]-1;
                if (0 <= new_coordinate and new_coordinate <= upper_limit) {
                    // We got a neighbor cell within the grid's limit
                    auto new_coordinates = coordinates;
                    new_coordinates[axis] += direction;
                    Cell & c = p_cells[p_grid->cell_index_at(new_coordinates)];
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
                        Cell * c = children.top();
                        children.pop();
                        GridCoordinates coordinates = subcell_coordinates[c->index];
                        if (coordinates[axis] == search_axis) {
                            if (c->data) {
                                neighbors.emplace_back(c);
                            } else {
                                for (auto & child : *c->childs) {
                                    children.emplace(&child);
                                }
                            }
                        }
                    }
                } else {
                    // No further childs, let's add the cell as a neighbor
                    neighbors.emplace_back(p);
                }
            }
        }
    }

    return neighbors;
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::populate_drawing_vectors()
{
    BEGIN_CLOCK;
    p_drawing_nodes_vector.resize(p_grid->number_of_nodes());
    p_drawing_edges_vector.resize(p_grid->number_of_edges()*2);
    // Reset the vector in case we got multiple initializations
    p_drawing_subdivided_edges_vector.resize(0);
    p_drawing_cells_vector.resize(0);
    p_drawing_subdivided_edges_vector.resize(p_regions.size());
    p_drawing_cells_vector.resize(p_regions.size());

    TICK;
    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_nodes(); ++i) {
        const auto p = p_grid->node(i);
        p_drawing_nodes_vector[i] = {p[0], p[1], p[2]};
    }

    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_edges(); ++i) {
        const auto edge = p_grid->edge(i);
        p_drawing_edges_vector[i*2] = {edge.node(0)[0], edge.node(0)[1], edge.node(0)[2]};
        p_drawing_edges_vector[i*2+1] = {edge.node(1)[0], edge.node(1)[1], edge.node(1)[2]};
    }

    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_cells(); ++i) {
        std::queue<std::pair<const Cell *, CellElement >> cells;
        cells.emplace(&p_cells[i], p_grid->cell_at(i));

        while (not cells.empty()) {
            const Cell * c = std::get<0>(cells.front());
            const CellElement & e = std::get<1>(cells.front());

            if (c->childs) {
                const auto & childs_elements = get_subcells(e);
                for (UNSIGNED_INTEGER_TYPE j = 0; j < (*c->childs).size(); ++j) {
                    cells.emplace(&(*c->childs)[j], childs_elements[j]);
                }
            } else {
                const auto & region_id = c->data->region_id;
                for (UNSIGNED_INTEGER_TYPE node_id = 0; node_id < ((unsigned) 1 << Dimension); ++node_id) {
                    p_drawing_cells_vector[region_id].emplace_back(e.node(node_id)[0], e.node(node_id)[1],
                                                                   e.node(node_id)[2]);
                }
                for (const auto &edge : CellElement::edges) {
                    p_drawing_subdivided_edges_vector[region_id].emplace_back(e.node(edge[0])[0], e.node(edge[0])[1],
                                                                   e.node(edge[0])[2]);
                    p_drawing_subdivided_edges_vector[region_id].emplace_back(e.node(edge[1])[0], e.node(edge[1])[1],
                                                                   e.node(edge[1])[2]);
                }
            }

            cells.pop();
        }
    }
    msg_info() << "Populating the drawing vectors in " << std::setprecision(3) << std::fixed
               << TOCK / 1000. / 1000.
               << " [ms]";
}


template <typename DataTypes>
std::pair<typename FictitiousGrid<DataTypes>::Coord, typename FictitiousGrid<DataTypes>::Coord>
FictitiousGrid<DataTypes>::compute_bbox_from(const SofaVecCoord & positions)
{
    Coord min = positions[0];
    Coord max = min;

    for (const Coord & p : positions) {
        for (size_t i = 0; i < Dimension; ++i) {
            if (p[i] < min[i])
                min[i] = p[i];
            if (p[i] > max[i])
                max[i] = p[i];
        }
    }

    return {min, max};
}

} // namespace SofaCaribou::GraphComponents::topology

#endif //SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL
