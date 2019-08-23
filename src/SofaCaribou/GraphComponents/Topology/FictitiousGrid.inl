#ifndef SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL
#define SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL

#include <SofaCaribou/GraphComponents/Topology/FictitiousGrid.h>

#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseTopology/MeshTopology.h>
#include <sofa/simulation/Node.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <chrono>

namespace SofaCaribou::GraphComponents::topology {

template <typename DataTypes>
FictitiousGrid<DataTypes>::FictitiousGrid()
        : d_n(initData(&d_n, SofaVecInt(),"n","Grid resolution [nx, ny, nz] where nx, ny or nz are the number of nodes in the x,y and z directions."))
        , d_min(initData(&d_min, SofaVecFloat(),"min","First corner node position of the grid's bounding box."))
        , d_max(initData(&d_max, SofaVecFloat(),"max","Second corner node position of the grid's bounding box."))
        , d_number_of_subdivision(initData(&d_number_of_subdivision,
                "maximum_number_of_subdivision_levels",
                0,
                "Number of subdivision levels of the boundary cells (one level split the cell in 8 subcells)."))
        , d_use_implicit_surface(initData(&d_use_implicit_surface,
                bool(false),
                "use_implicit_surface",
                "Use an implicit surface instead of a tessellated surface. If true, the callback function is_inside must be defined."))
        , d_surface_positions(initData(&d_surface_positions, SofaVecCoord(),
                "surface_positions",
                "Position vector of nodes contained in the surface boundary of the immersed object."))
        , d_surface_edges(initData(&d_surface_edges,
                "surface_edges",
                "List of edges (ex: [e1p1 e1p2 e2p1 e2p2 ...])."))
        , d_surface_triangles(initData(&d_surface_triangles,
                "surface_triangles",
                "List of triangles (ex: [t1p1 t1p2 t1p3 t2p1 t2p2 t2p3 ...])."))
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
                        sofa::helper::WriteAccessor<Data<sofa::helper::vector<Edge>>> w_edges = d_surface_edges;
                        for (unsigned int i = 0; i < topology->getNbEdges(); ++i) {
                            w_edges.push_back(topology->getEdge(i));
                        }
                        msg_info() << "Automatically found " << topology->getNbEdges() << " edges in the container '"
                                   << topology->getPathName() << "'.";
                    } else {
                        sofa::helper::WriteAccessor<Data<sofa::helper::vector<Triangle >>> w_triangles = d_surface_triangles;
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
                        sofa::helper::WriteAccessor<Data<sofa::helper::vector<Edge>>> w_edges = d_surface_edges;
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
                        sofa::helper::WriteAccessor<Data<sofa::helper::vector<Triangle >>> w_triangles = d_surface_triangles;
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

    p_node_types.resize(p_grid->number_of_nodes(), Type::Undefined);
    p_cells_types.resize(p_grid->number_of_cells(), Type::Undefined);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    if (d_use_implicit_surface.getValue() and p_implicit_test_callback) {
        compute_cell_types_from_implicit_surface();
    } else {
        compute_cell_types_from_explicit_surface();
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    msg_info() << "Found the boundary cells in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]";
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

template <typename DataTypes>
void FictitiousGrid<DataTypes>::draw(const sofa::core::visual::VisualParams* vparams)
{
    using Color = sofa::defaulttype::Vec4f;
    if (!p_grid)
        return;

    std::vector<sofa::defaulttype::Vector3> points_inside;
    std::vector<sofa::defaulttype::Vector3> points_boundary;
    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_nodes(); ++i) {
        const auto p = p_grid->node(i);
        if (p_node_types[i] == Type::Inside)
            points_inside.push_back({p[0], p[1], p[2]});
        else if (p_node_types[i] == Type::Boundary)
            points_boundary.push_back({p[0], p[1], p[2]});
    }
    vparams->drawTool()->drawPoints(points_inside, 5, Color(0, 1, 0, 1));
    vparams->drawTool()->drawPoints(points_boundary, 5, Color(1, 0, 0, 1));


    std::vector<sofa::defaulttype::Vector3> edge_points(p_grid->number_of_edges()*2);
    for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_edges(); ++i) {
        const auto edge = p_grid->edge(i);
        edge_points[i*2] = {edge.node(0)[0], edge.node(0)[1], edge.node(0)[2]};
        edge_points[i*2+1] = {edge.node(1)[0], edge.node(1)[1], edge.node(1)[2]};
    }

    vparams->drawTool()->drawLines(edge_points, 0.5, Color(1,1,1, 1));

}

} // namespace SofaCaribou::GraphComponents::topology

#endif //SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL
