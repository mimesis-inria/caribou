#ifndef SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL
#define SOFACARIBOU_GRAPHCOMPONENTS_TOPOLOGY_FICTITIOUSGRID_INL

#include <SofaCaribou/GraphComponents/Topology/FictitiousGrid.h>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/Node.h>
#include <sofa/core/behavior/MechanicalState.h>

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
                bool(true),
                "use_implicit_surface",
                "Use an implicit surface instead of a tessellated surface. If true, the callback function is_inside must be defined."))
        , d_surface_positions(initData(&d_surface_positions, SofaVecCoord(),
                "surface_positions",
                "Position vector of nodes contained in the surface boundary of the immersed object."))
        , d_surface_topology_container(initLink(
                "surface_topology_container",
                "Topology that contains the elements on the surface boundary of the immersed object."))
{
}

template <typename DataTypes>
void FictitiousGrid<DataTypes>::init() {
    using namespace sofa::core::topology;
    using namespace sofa::core::behavior;

    if (not d_use_implicit_surface.getValue()) {
        // Fetch the surface topology container
        if (not d_surface_topology_container.get()) {
            auto containers = this->getContext()->template getObjects<TopologyContainer>(BaseContext::Local);
            auto node = static_cast<const sofa::simulation::Node *> (this->getContext());
            if (containers.empty() || containers.size() > 1) {

                std::vector<TopologyContainer *> surface_containers;
                for (TopologyContainer *c : containers) {
                    if (Dimension == 2 and c->getTopologyType() == sofa::core::topology::EDGE)
                        surface_containers.emplace_back(c);
                    else if (Dimension == 3 and (c->getTopologyType() == sofa::core::topology::TRIANGLE ||
                                                 c->getTopologyType() == sofa::core::topology::QUAD))
                        surface_containers.emplace_back(c);
                }

                if (surface_containers.empty()) {
                    msg_error()
                        << "No topology container were found in the context node '" << node->getPathName() << "'. "
                        << "Please set the parameter '" << d_surface_topology_container.getName()
                        << "' to the path of the topology container that "
                        << "contains the surface boundary mesh of the immersed object, or place a topology container in the current context. "
                        << "Note that in two dimensions, the surface topology must be an edge topology container. In three dimensions, "
                        << "the container must be either a triangle or a quad topology container.";
                    return;
                } else if (surface_containers.size() > 1) {
                    msg_error()
                        << "Multiple topology containers were found in the context node '" << node->getPathName()
                        << "'."
                        << "Please set the parameter '" << d_surface_topology_container.getName()
                        << "' to the path of the topology container that "
                        << "contains the surface boundary mesh of the immersed object. "
                        << "Note that in two dimensions, the surface topology must be an edge topology container. In three dimensions, "
                        << "the container must be either a triangle or a quad topology container.";
                    return;
                } else {
                    d_surface_topology_container.set(surface_containers[0]);
                    msg_info() << "Automatically found the topology container at '"
                               << surface_containers[0]->getPathName() << "'.";
                }
            } else {
                d_surface_topology_container.set(containers[0]);
                msg_info() << "Automatically found the topology container at '" << containers[0]->getPathName() << "'.";
            }
        }

        TopologyContainer *topology = d_surface_topology_container.get();
        if (Dimension == 2 and topology->getNbEdges() == 0) {
            msg_warning() << "The topology container '" << topology->getPathName()
                          << "' does not contain any edge elements.";
        } else if (Dimension == 3 and (topology->getNbTriangles() == 0 and topology->getNbQuads() == 0)) {
            msg_warning() << "The topology container '" << topology->getPathName()
                          << "' does not contain any triangle or quad elements.";
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
                    auto surface_context = topology->getContext();
                    auto surface_node = dynamic_cast<const sofa::simulation::Node *> (surface_context);
                    auto mechanical_states = surface_context->template getObjects<MechanicalState<DataTypes>>(
                        BaseContext::Local);
                    if (mechanical_states.empty()) {
                        msg_error() << "No positions were set in the '" << d_surface_positions.getName()
                                    << "' data parameter, "
                                    << "and no mechanical objects were found in the surface topology context node ('"
                                    << surface_node->getPathName() << "').";
                        return;
                    } else if (mechanical_states.size() > 1) {
                        msg_error() << "No positions were set in the '" << d_surface_positions.getName()
                                    << "' data parameter, "
                                    << "and more than one mechanical objects were found in the surface topology context node ('"
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
                                msg_error() << "A mechanical state was automatically found at '" << mo->getPathName()
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
    compute_cell_types();
}

template <typename DataTypes>
void
FictitiousGrid<DataTypes>::compute_cell_types()
{
    if (!p_grid or p_grid->number_of_nodes() == 0)
        return;

    p_node_type.resize(p_grid->number_of_nodes(), Type::Undefined);

    if (d_use_implicit_surface.getValue() and p_implicit_test_callback) {
        for (UNSIGNED_INTEGER_TYPE i = 0; i < p_grid->number_of_nodes(); ++i) {
            const auto t = p_implicit_test_callback(p_grid->node(i));
            if (t < 0)
                p_node_type[i] = Type::Inside;
            else if (t > 0)
                p_node_type[i] = Type::Outside;
            else
                p_node_type[i] = Type::Boundary;
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
        if (p_node_type[i] == Type::Inside)
            points_inside.push_back({p[0], p[1], p[2]});
        else if (p_node_type[i] == Type::Boundary)
            points_boundary.push_back({p[0], p[1], p[2]});
    }
    vparams->drawTool()->drawPoints(points_inside, 1, Color(0, 1, 0, 1));
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
