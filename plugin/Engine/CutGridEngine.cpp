#undef NDEBUG
#include "CutGridEngine.h"
#include "../Helper/Hexahedron.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/visual/DrawToolGL.h>

namespace sofa
{

namespace caribou
{

namespace engine
{

CutGridEngine::CutGridEngine()
    : d_grid_topology(initLink("grid_topology", "[IN] Grid topology that contains the regular hexahedrons"))
    , d_surface_topology(initLink("surface_topology", "[IN] Mesh topology that contains the surface elements"))
    , d_showHexahedrons(initData(&d_showHexahedrons, "showHexahedrons","[IN] \"None\", \"Inside\", \"Outside\", \"Boundary\" or \"All\""))
    , d_showTriangles(initData(&d_showTriangles, false, "showTriangles", "[IN] Show the triangles created in hexahedrons that lie on the boundary."))
    , d_showInsideColor(initData(&d_showInsideColor, Color(1, 0, 0, 1), "showInsideColor", "[IN] Color of the hexahedrons that lie inside the surface"))
    , d_showOutsideColor(initData(&d_showOutsideColor, Color(0, 1, 0, 1), "showOutsideColor", "[IN] Color of the hexahedrons that lie outside the surface"))
    , d_showBoundaryColor(initData(&d_showBoundaryColor, Color(0, 0, 1, 1), "showBoundaryColor", "[IN] Color of the hexahedrons that lie on the boundary of the surface"))
    , d_showTrianglesColor(initData(&d_showTrianglesColor, Color(1, 0, 1, 1), "showTrianglesColor", "[IN] Color of the triangles created in hexahedrons that lie on the boundary."))
    , d_points_flags(initData(&d_points_flags, "points_flags", "[OUT] Flags (Outside -1, Boundary 0 or Inside 1) of the points"))
    , d_hexahedrons_flags(initData(&d_hexahedrons_flags, "hexahedrons_flags", "[OUT] Flags (Outside -1, Boundary 0 or Inside 1) of the hexahedrons"))
    , d_triangle_positions(initData(&d_triangle_positions, "triangle_positions", "[OUT] List of triangles vertices positions."))
    , d_triangles(initData(&d_triangles, "triangles", "[OUT] List of created triangles per hexahedron.", false))
{
    d_points_flags.setReadOnly(true);
    d_hexahedrons_flags.setReadOnly(true);
    d_triangle_positions.setReadOnly(true);
    d_triangles.setReadOnly(true);

    sofa::helper::WriteAccessor<Data<sofa::helper::OptionsGroup>> showHexahedronsOptions = d_showHexahedrons;
    showHexahedronsOptions->setNames(5,"NONE", "INSIDE", "OUTSIDE", "BOUNDARY", "ALL");
    showHexahedronsOptions->setSelectedItem(0);
}

void CutGridEngine::init()
{
    Inherit::init();

    // Inputs
    addInput(&d_showHexahedrons);
    addInput(&d_showTriangles);
    addInput(&d_showInsideColor);
    addInput(&d_showOutsideColor);
    addInput(&d_showBoundaryColor);
    addInput(&d_showTrianglesColor);

    // Outputs
    addOutput(&d_points_flags);
    addOutput(&d_hexahedrons_flags);
    addOutput(&d_triangle_positions);
    addOutput(&d_triangles);

    setDirtyValue();
}

void CutGridEngine::update()
{
    const auto number_of_hexahedrons = d_grid_topology->getNbHexahedra();
    const auto number_of_points = d_grid_topology->getNbPoints();

    using EdgeId = sofa::caribou::helper::hexahedron::EdgeId;
    using AlphaValue = sofa::caribou::helper::hexahedron::AlphaValue;

    cleanDirty();

    sofa::helper::WriteOnlyAccessor<Data<sofa::helper::vector<Flag>>> points_flags = d_points_flags;
    sofa::helper::WriteOnlyAccessor<Data<sofa::helper::vector<Flag>>> hexahedrons_flags = d_hexahedrons_flags;
    sofa::helper::WriteOnlyAccessor<Data<sofa::helper::vector<Coord>>> triangle_positions = d_triangle_positions;
    sofa::helper::WriteOnlyAccessor<Data<sofa::helper::vector<sofa::helper::vector<Triangle>>>> triangles = d_triangles;

    points_flags.resize(number_of_points);
    hexahedrons_flags.resize(number_of_hexahedrons);
    triangles.resize(number_of_hexahedrons);

    for (size_t i = 0; i < points_flags.size(); ++i) {
        const Coord & p = d_grid_topology->getPoint(i);
        points_flags[i] = getFlag(p);
    }


    std::map< Coord, PointID> vertices;
    size_t nb_non_boundary_with_triangles = 0;
    for (size_t i = 0; i < hexahedrons_flags.size(); ++i) {
        const auto & hexa = d_grid_topology->getHexa(i);
        triangles[i].clear();

        std::array<Coord, 8> nodes;
        for (unsigned char j = 0; j < 8; ++j)
            nodes[j] = d_grid_topology->getPoint(hexa[j]);

        sofa::helper::vector<std::array<Coord, 3>> cube_triangles = sofa::caribou::helper::hexahedron::triangulate_interior(
                nodes,
                [this, nodes](const unsigned char & node_id) {
                    Flag flag = getFlag(nodes[node_id]);
                    return flag == Flag::Inside || flag == Flag::Boundary;
                },
                [nodes](const EdgeId & edge_id) {
                    constexpr Real x0 = 0, x1 = 10, radius = 0.5;
                    Coord p0 = nodes[sofa::caribou::helper::hexahedron::edges[edge_id][0]];
                    Coord p1 = nodes[sofa::caribou::helper::hexahedron::edges[edge_id][1]];

                    if (p0[0] < x0)
                        p0[0] = x0;
                    else if (p0[0] > x1)
                        p0[0] = x1;

                    if (p1[0] < x0)
                        p1[0] = x0;
                    else if (p1[0] > x1)
                        p1[0] = x1;

                    auto isovalue_0 = (AlphaValue) (radius - sqrt(p0[1]*p0[1] + p0[2]*p0[2]));
                    auto isovalue_1 = (AlphaValue) (radius - sqrt(p1[1]*p1[1] + p1[2]*p1[2]));

                    AlphaValue alpha = - isovalue_0 / (isovalue_1 - isovalue_0);
                    return alpha;
                });

        for (const auto & t : cube_triangles) {
            Triangle triangle;
            for (unsigned char j = 0; j < 3; ++j) {
                auto iter = vertices.find(t[j]);
                if (iter != vertices.end()) {
                    triangle[j] = iter->second;
                } else {
                    PointID id = vertices.size();
                    vertices.insert(std::make_pair(t[j], id));
                    triangle[j] = id;
                }
            }
            triangles[i].push_back(triangle);
        }

        hexahedrons_flags[i] = getFlag(hexa);

        if (hexahedrons_flags[i] == Flag::Boundary && cube_triangles.empty()) {
            hexahedrons_flags[i] = Flag::Outside;
        } else if (hexahedrons_flags[i] != Flag::Boundary && !cube_triangles.empty()) {
            nb_non_boundary_with_triangles++;
        }
    }

    if (nb_non_boundary_with_triangles > 0) {
        msg_error() << std::to_string(nb_non_boundary_with_triangles) << " hexahedrons of type inside or outside but with some triangles.";
    }

    triangle_positions.resize(vertices.size());
    for (const auto & v : vertices) {
        triangle_positions[v.second] = v.first;
    }

}

CutGridEngine::Flag CutGridEngine::getFlag(const Coord & p) const {

    constexpr Real x0 = 0, x1 = 10, radius = 0.5;

    const Real & x = p[0];
    const Real & y = p[1];
    const Real & z = p[2];

    if (x < x0 || x > x1)
        return Flag::Outside;

    Real dist = sqrt(y*y + z*z);
    if (dist > radius)
        return Flag::Outside;

//    if (abs(dist - radius) < epsilon && ((abs(x - x0) < epsilon) || (abs(x - x1) < epsilon)))
//        return Flag::Boundary;

    return Flag::Inside;
}

CutGridEngine::Flag CutGridEngine::getFlag(const Hexa & hexa) const {
    bool all_inside = true, all_outside = true;
    for (unsigned char j = 0; j < 8; ++j) {
        Flag flag = getFlag(d_grid_topology->getPoint(hexa[j]));
        switch (flag) {
            case Flag::Outside: all_inside = false; break;
            case Flag::Inside:  all_outside = false; break;
            case Flag::Boundary:
            default:
                all_inside = false;
                all_outside = false;
        }
    }

    if (all_inside)
        return Flag::Inside;
    else if(all_outside)
        return Flag::Outside;
    else
        return Flag::Boundary;
}

void CutGridEngine::draw(const core::visual::VisualParams* vparams)
{
    using Vector3 = core::visual::DrawTool::Vector3;
    const auto drawTool = static_cast<sofa::core::visual::DrawToolGL *> (vparams->drawTool());

    unsigned int showHexahedron = d_showHexahedrons.getValue().getSelectedId();
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Flag>>> hexahedrons_flags = d_hexahedrons_flags;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Coord>>> triangle_positions = d_triangle_positions;
    sofa::helper::ReadAccessor<Data<sofa::helper::vector<sofa::helper::vector<Triangle>>>> triangles = d_triangles;

    if (showHexahedron > 0 /*NONE*/) {
        std::vector<Vector3> points_inside;
        std::vector<Vector3> edges_inside;
        std::vector<Vector3> points_outside;
        std::vector<Vector3> edges_outside;
        std::vector<Vector3> points_boundary;
        std::vector<Vector3> edges_boundary;
        const auto number_of_hexahedrons = d_grid_topology->getNbHexahedra();
        for (size_t i = 0; i < number_of_hexahedrons; ++i) {
            if (showHexahedron == 4 /*ALL*/
            || (showHexahedron == 1 /*INSIDE*/ && hexahedrons_flags[i] == Flag::Inside)
            || (showHexahedron == 2 /*OUTSIDE*/ && hexahedrons_flags[i] == Flag::Outside)
            || (showHexahedron == 3 /*BOUNDARY*/ && hexahedrons_flags[i] == Flag::Boundary)) {
                for (unsigned char j = 0; j < 8; ++j) {
                    Coord point = d_grid_topology->getPoint(d_grid_topology->getHexa(i)[j]);
                    switch (hexahedrons_flags[i]) {
                        case Flag::Inside: points_inside.push_back(point); break;
                        case Flag::Outside: points_outside.push_back(point); break;
                        case Flag::Boundary: points_boundary.push_back(point); break;
                        default: break;
                    }
                }

                for (const auto & edge : sofa::caribou::helper::hexahedron::edges) {
                    Coord point1 = d_grid_topology->getPoint(d_grid_topology->getHexa(i)[edge[0]]);
                    Coord point2 = d_grid_topology->getPoint(d_grid_topology->getHexa(i)[edge[1]]);
                    switch (hexahedrons_flags[i]) {
                        case Flag::Inside: edges_inside.push_back(point1); edges_inside.push_back(point2); break;
                        case Flag::Outside: edges_outside.push_back(point1); edges_outside.push_back(point2); break;
                        case Flag::Boundary: edges_boundary.push_back(point1); edges_boundary.push_back(point2); break;
                        default: break;
                    }
                }
            }
        }

        Color insideColor = d_showInsideColor.getValue();
        Color outsideColor = d_showOutsideColor.getValue();
        Color boundaryColor = d_showBoundaryColor.getValue();

        drawTool->drawHexahedra(points_inside, insideColor);
        drawTool->drawHexahedra(points_outside, outsideColor);
        drawTool->drawHexahedra(points_boundary, boundaryColor);

        insideColor[3] = 1;
        outsideColor[3] = 1;
        boundaryColor[3] = 1;

        drawTool->drawLines(edges_inside, 1, insideColor);
        drawTool->drawLines(edges_outside, 1, outsideColor);
        drawTool->drawLines(edges_boundary, 1, boundaryColor);
    }

    if (d_showTriangles.getValue()) {
        std::vector<Vector3> points;
        std::vector<Vector3> edges_points;
        const auto number_of_hexahedrons = d_grid_topology->getNbHexahedra();
        for (size_t i = 0; i < number_of_hexahedrons; ++i) {
            for (const auto &triangle : triangles[i]) {
                points.push_back(triangle_positions[triangle[0]]);
                edges_points.push_back(triangle_positions[triangle[0]]);

                points.push_back(triangle_positions[triangle[1]]);
                edges_points.push_back(triangle_positions[triangle[1]]);
                edges_points.push_back(triangle_positions[triangle[1]]);

                points.push_back(triangle_positions[triangle[2]]);
                edges_points.push_back(triangle_positions[triangle[2]]);
                edges_points.push_back(triangle_positions[triangle[2]]);

                edges_points.push_back(triangle_positions[triangle[0]]);
            }
        }
        drawTool->drawTriangles(points, d_showTrianglesColor.getValue());
        drawTool->drawLines(edges_points, 1, d_showTrianglesColor.getValue());
    }

//    sofa::helper::ReadAccessor<Data<sofa::helper::vector<Flag>>> points_flags = d_points_flags;
//    sofa::helper::ReadAccessor<Data<VecCoord>> positions = this->getMState()->readPositions();
//
//    if (points_flags.size() != positions.size())
//        return;
//
//    sofa::helper::vector<Vector3> points (positions.size());
//    sofa::helper::vector<Color> colors (positions.size());
//
//    for (size_t i = 0; i < points_flags.size(); ++i) {
//        const Flag & flag = points_flags[i];
//        switch (flag) {
//            case Flag::Outside:
//                colors[i] = Color(1, 0, 0, 1);
//                break;
//            case Flag::Boundary:
//                colors[i] = Color(0, 1, 0, 1);
//                break;
//            case Flag::Inside:
//                colors[i] = Color(1, 1, 1, 1);
//                break;
//            default:
//                colors[i] = Color(0, 0, 0, 0);
//        }
//
//        points[i] = positions[i];
//    }
//
//    vparams->drawTool()->drawPoints(points, 5, colors);
}

static int CutGridEngineClass = core::RegisterObject("Caribou cut grid engine")
        .add< CutGridEngine >(true)
;

} // namespace topology

} // namespace caribou

} // namespace sofa