#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/testing/BaseTest.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/simulation/graph/DAGSimulation.h>
#include <sofa/simulation/graph/SimpleApi.h>
#include <sofa/helper/system/PluginManager.h>
#include <sofa/simulation/Node.h>
DISABLE_ALL_WARNINGS_BEGIN

#include <SofaCaribou/Topology/CaribouTopology[Quad].h>
#include <SofaCaribou/Topology/CaribouTopology[Quad8].h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron].h>
#include <SofaCaribou/Topology/CaribouTopology[Tetrahedron10].h>
#include <SofaCaribou/Topology/CaribouTopology[Triangle].h>
#include <SofaCaribou/Topology/CaribouTopology[Triangle6].h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron20].h>
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/IO/VTKReader.h>

#include "../sofacaribou_test.h"

using sofa::helper::system::PluginManager ;
using namespace sofa::simulation;
using namespace sofa::simpleapi;
using namespace sofa::helper::logging;
using namespace sofa::testing;

TEST(CaribouTopology, QuadLinear2DAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_2D, PointID>::MeshType;
    using Domain = Mesh::Domain<Quad<_2D>, PointID>;
    auto reader = io::VTKReader<_2D, PointID>::Read(executable_directory_path + "/meshes/2D_quad_linear.vtk");

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 25);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the quad surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Quad<_2D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Quad_2D"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 4>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 16); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 4; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, QuadLinear2DFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_2D, PointID>::MeshType;
    auto reader = io::VTKReader<_2D, PointID>::Read(executable_directory_path + "/meshes/2D_quad_linear.vtk");
    using Domain = Mesh::Domain<Quad<_2D>, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 25);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the quad surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec2Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}, {"template", "Vec2d"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Quad<_2D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Quad_2D"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 4>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 4; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * quad_domain = topo->domain();
    ASSERT_NE(quad_domain, nullptr);

    FLOATING_POINT_TYPE area = 0;
    for (UNSIGNED_INTEGER_TYPE quad_id = 0; quad_id < quad_domain->number_of_elements(); ++quad_id) {
        auto quad = quad_domain->element(quad_id);
        for (const auto & g : quad.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto detJ = quad.jacobian(x).determinant();
            area += w*abs(detJ);
        }
    }
    EXPECT_NEAR(area, 100., 1e-3);
}


TEST(CaribouTopology, QuadQuadratic2DAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_2D, PointID>::MeshType;
    using Domain = Mesh::Domain<Quad8<_2D>, PointID>;
    auto reader = io::VTKReader<_2D, PointID>::Read(executable_directory_path + "/meshes/2D_quad_quadratic.vtk");

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 65);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the quad surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Quad8<_2D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Quad8_2D"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 8>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 16); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 8; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, QuadQuadratic2DFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_2D, PointID>::MeshType;
    auto reader = io::VTKReader<_2D, PointID>::Read(executable_directory_path + "/meshes/2D_quad_quadratic.vtk");
    using Domain = Mesh::Domain<Quad8<_2D>, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 65);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the quad surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec2Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}, {"template", "Vec2d"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Quad8<_2D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Quad8_2D"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 8>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 8; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * quad_domain = topo->domain();
    ASSERT_NE(quad_domain, nullptr);

    FLOATING_POINT_TYPE area = 0;
    for (UNSIGNED_INTEGER_TYPE quad_id = 0; quad_id < quad_domain->number_of_elements(); ++quad_id) {
        auto quad = quad_domain->element(quad_id);
        for (const auto & g : quad.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto detJ = quad.jacobian(x).determinant();
            area += w*abs(detJ);
        }
    }
    EXPECT_NEAR(area, 100., 1e-3);
}


TEST(CaribouTopology, QuadLinear3DAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    using Domain = Mesh::Domain<Quad<_3D>, PointID>;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_quad_linear.vtk");

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 25);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the quad surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Quad<_3D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Quad"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 4>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 16); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 4; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, QuadLinear3DFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_quad_linear.vtk");
    using Domain = Mesh::Domain<Quad<_3D>, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 25);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the quad surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}, {"template", "Vec3d"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1], p[2]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Quad<_3D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Quad"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 4>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 4; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * quad_domain = topo->domain();
    ASSERT_NE(quad_domain, nullptr);

    FLOATING_POINT_TYPE area = 0;
    for (UNSIGNED_INTEGER_TYPE quad_id = 0; quad_id < quad_domain->number_of_elements(); ++quad_id) {
        auto quad = quad_domain->element(quad_id);
        for (const auto & g : quad.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto J = quad.jacobian(x);
            const auto detJ = J.col(0).cross(J.col(1)).norm();
            area += w*abs(detJ);
        }
    }
    EXPECT_NEAR(area, 100., 1e-3);
}


TEST(CaribouTopology, QuadQuadratic3DAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    using Domain = Mesh::Domain<Quad8<_3D>, PointID>;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_quad_quadratic.vtk");

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 65);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the quad surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Quad8<_3D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Quad8"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 8>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 16); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 8; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, QuadQuadratic3DFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_quad_quadratic.vtk");
    using Domain = Mesh::Domain<Quad8<_3D>, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 65);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the quad surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}, {"template", "Vec3d"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1], p[2]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Quad8<_3D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Quad8"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 8>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 8; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * quad_domain = topo->domain();
    ASSERT_NE(quad_domain, nullptr);

    FLOATING_POINT_TYPE area = 0;
    for (UNSIGNED_INTEGER_TYPE quad_id = 0; quad_id < quad_domain->number_of_elements(); ++quad_id) {
        auto quad = quad_domain->element(quad_id);
        for (const auto & g : quad.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto J = quad.jacobian(x);
            const auto detJ = J.col(0).cross(J.col(1)).norm();
            area += w*abs(detJ);
        }
    }
    EXPECT_NEAR(area, 100., 1e-3);
}


TEST(CaribouTopology, TriangleLinear2DAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_2D, PointID>::MeshType;
    using Domain = Mesh::Domain<Triangle<_2D>, PointID>;
    auto reader = io::VTKReader<_2D, PointID>::Read(executable_directory_path + "/meshes/2D_triangle_linear.vtk");

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 25);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the triangle surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Triangle<_2D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Triangle_2D"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 3>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 32); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 3; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, TriangleLinear2DFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_2D, PointID>::MeshType;
    auto reader = io::VTKReader<_2D, PointID>::Read(executable_directory_path + "/meshes/2D_triangle_linear.vtk");
    using Domain = Mesh::Domain<Triangle<_2D>, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 25);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the triangle surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec2Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}, {"template", "Vec2d"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Triangle<_2D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Triangle_2D"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 3>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 3; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * triangle_domain = topo->domain();
    ASSERT_NE(triangle_domain, nullptr);

    FLOATING_POINT_TYPE area = 0;
    for (UNSIGNED_INTEGER_TYPE triangle_id = 0; triangle_id < triangle_domain->number_of_elements(); ++triangle_id) {
        auto triangle = triangle_domain->element(triangle_id);
        for (const auto & g : triangle.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto detJ = triangle.jacobian(x).determinant();
            area += w*abs(detJ);
        }
    }
    EXPECT_NEAR(area, 100., 1e-3);
}


TEST(CaribouTopology, TriangleQuadratic2DAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_2D, PointID>::MeshType;
    using Domain = Mesh::Domain<Triangle6<_2D>, PointID>;
    auto reader = io::VTKReader<_2D, PointID>::Read(executable_directory_path + "/meshes/2D_triangle_quadratic.vtk");

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 81);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the triangle surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Triangle6<_2D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Triangle6_2D"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 6>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 32); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 6; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, TriangleQuadratic2DFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_2D, PointID>::MeshType;
    auto reader = io::VTKReader<_2D, PointID>::Read(executable_directory_path + "/meshes/2D_triangle_quadratic.vtk");
    using Domain = Mesh::Domain<Triangle6<_2D>, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 81);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the triangle surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec2Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}, {"template", "Vec2d"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Triangle6<_2D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Triangle6_2D"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 6>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 6; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * triangle_domain = topo->domain();
    ASSERT_NE(triangle_domain, nullptr);

    FLOATING_POINT_TYPE area = 0;
    for (UNSIGNED_INTEGER_TYPE triangle_id = 0; triangle_id < triangle_domain->number_of_elements(); ++triangle_id) {
        auto triangle = triangle_domain->element(triangle_id);
        for (const auto & g : triangle.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto detJ = triangle.jacobian(x).determinant();
            area += w*abs(detJ);
        }
    }
    EXPECT_NEAR(area, 100., 1e-3);
}


TEST(CaribouTopology, TriangleLinear3DAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    using Domain = Mesh::Domain<Triangle<_3D>, PointID>;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_triangle_linear.vtk");

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 25);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the triangle surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Triangle<_3D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Triangle"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 3>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 32); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 3; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, TriangleLinear3DFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_triangle_linear.vtk");
    using Domain = Mesh::Domain<Triangle<_3D>, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 25);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the triangle surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1], p[2]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Triangle<_3D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Triangle"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 3>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 3; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * triangle_domain = topo->domain();
    ASSERT_NE(triangle_domain, nullptr);

    FLOATING_POINT_TYPE area = 0;
    for (UNSIGNED_INTEGER_TYPE triangle_id = 0; triangle_id < triangle_domain->number_of_elements(); ++triangle_id) {
        auto triangle = triangle_domain->element(triangle_id);
        for (const auto & g : triangle.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto J = triangle.jacobian(x);
            const auto detJ = J.col(0).cross(J.col(1)).norm();
            area += w*detJ;
        }
    }
    EXPECT_NEAR(area, 100., 1e-3);
}


TEST(CaribouTopology, TriangleQuadratic3DAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    using Domain = Mesh::Domain<Triangle6<_3D>, PointID>;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_triangle_quadratic.vtk");

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 81);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the triangle surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Triangle6<_3D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Triangle6"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 6>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 32); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 6; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, TriangleQuadratic3DFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_triangle_quadratic.vtk");
    using Domain = Mesh::Domain<Triangle6<_3D>, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 81);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the segment contour, second is the triangle surface domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1], p[2]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Triangle6<_3D>> *> (
            createObject(root, "CaribouTopology", {{"template", "Triangle6"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 6>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 6; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * triangle_domain = topo->domain();
    ASSERT_NE(triangle_domain, nullptr);

    FLOATING_POINT_TYPE area = 0;
    for (UNSIGNED_INTEGER_TYPE triangle_id = 0; triangle_id < triangle_domain->number_of_elements(); ++triangle_id) {
        auto triangle = triangle_domain->element(triangle_id);
        for (const auto & g : triangle.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto J = triangle.jacobian(x);
            const auto detJ = J.col(0).cross(J.col(1)).norm();
            area += w*detJ;
        }
    }
    EXPECT_NEAR(area, 100., 1e-3);
}


TEST(CaribouTopology, TetrahedronLinearAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_tetrahedron_linear.vtk");
    using Domain = Mesh::Domain<Tetrahedron, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 206);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the triangular surface, second is the tetra volumetric domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Tetrahedron> *> (
        createObject(root, "CaribouTopology", {{"template", "Tetrahedron"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 4>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 681); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 4; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
            << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
            << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, TetrahedronLinearFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_tetrahedron_linear.vtk");
    using Domain = Mesh::Domain<Tetrahedron, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 206);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the triangular surface, second is the tetra volumetric domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1], p[2]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Tetrahedron> *> (
            createObject(root, "CaribouTopology", {{"template", "Tetrahedron"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 4>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 4; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    double pi = 3.14159265358979323846;
    double r = 5;
    const auto * tetra_domain = topo->domain();
    ASSERT_NE(tetra_domain, nullptr);

    FLOATING_POINT_TYPE volume = 0;
    for (UNSIGNED_INTEGER_TYPE tetra_id = 0; tetra_id < tetra_domain->number_of_elements(); ++tetra_id) {
        auto tetra = tetra_domain->element(tetra_id);
        for (const auto & g : tetra.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto J = tetra.jacobian(x);
            const auto detJ = J.determinant();
            volume += w*abs(detJ);
        }
    }
    double exact_volume = 4/3. * pi * r*r*r;
    EXPECT_LE( abs((volume-exact_volume)/exact_volume), 0.05); // 5% error max
}


TEST(CaribouTopology, TetrahedronQuadraticAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_tetrahedron_quadratic.vtk");
    using Domain = Mesh::Domain<Tetrahedron10, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 1250);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the triangular surface, second is the tetra volumetric domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Tetrahedron10> *> (
            createObject(root, "CaribouTopology", {{"template", "Tetrahedron10"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 10>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 681); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 10; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, TetrahedronQuadraticFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_tetrahedron_quadratic.vtk");
    using Domain = Mesh::Domain<Tetrahedron10, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 1250);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the triangular surface, second is the tetra volumetric domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1], p[2]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Tetrahedron10> *> (
            createObject(root, "CaribouTopology", {{"template", "Tetrahedron10"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 10>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 10; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    double pi = 3.14159265358979323846;
    double r = 5;
    const auto * tetra_domain = topo->domain();
    ASSERT_NE(tetra_domain, nullptr);

    FLOATING_POINT_TYPE volume = 0;
    for (UNSIGNED_INTEGER_TYPE tetra_id = 0; tetra_id < tetra_domain->number_of_elements(); ++tetra_id) {
        auto tetra = tetra_domain->element(tetra_id);
        for (const auto & g : tetra.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto J = tetra.jacobian(x);
            const auto detJ = J.determinant();
            volume += w*abs(detJ);
        }
    }
    double exact_volume = 4/3. * pi * r*r*r;
    EXPECT_LE( abs((volume-exact_volume)/exact_volume), 0.05); // 5% error max
}


TEST(CaribouTopology, HexahedronLinearAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_hexahedron_linear.vtk");
    using Domain = Mesh::Domain<Hexahedron, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 125);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the quad surface, second is the hexa volumetric domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Hexahedron> *> (
            createObject(root, "CaribouTopology", {{"template", "Hexahedron"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 8>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 64); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 8; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, HexahedronLinearFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_hexahedron_linear.vtk");
    using Domain = Mesh::Domain<Hexahedron, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 125);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the quad surface, second is the hexa volumetric domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1], p[2]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Hexahedron> *> (
            createObject(root, "CaribouTopology", {{"template", "Hexahedron"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 8>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 8; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * hexa_domain = topo->domain();
    ASSERT_NE(hexa_domain, nullptr);

    FLOATING_POINT_TYPE volume = 0;
    for (UNSIGNED_INTEGER_TYPE hexa_id = 0; hexa_id < hexa_domain->number_of_elements(); ++hexa_id) {
        auto hexa = hexa_domain->element(hexa_id);
        for (const auto & g : hexa.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto J = hexa.jacobian(x);
            const auto detJ = J.determinant();
            volume += w*abs(detJ);
        }
    }
    double exact_volume = 10*10*10;
    EXPECT_LE( abs((volume-exact_volume)/exact_volume), 0.001); // 1% error max
}


TEST(CaribouTopology, HexahedronQuadraticAttachDomain) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_hexahedron_quadratic.vtk");
    using Domain = Mesh::Domain<Hexahedron20, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 425);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the quad surface, second is the hexa volumetric domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Hexahedron20> *> (
            createObject(root, "CaribouTopology", {{"template", "Hexahedron20"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Attach the domain
    topo->attachDomain(domain);

    // Make sure the data parameter `indices` has been filled-up correctly
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 20>>>;
    auto indices = ReadAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    EXPECT_EQ(indices.size(), 64); // Number of elements

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        const auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 20; ++node_id) {
            EXPECT_EQ(element_indices[node_id], domain->element_indices(element_id)[node_id])
                                << "The indices of the "<<node_id<<"th node of the element #"<<element_id<<" is "<< element_indices[node_id]
                                << " while it should be "<<domain->element_indices(element_id)[node_id];
        }
    }
}

TEST(CaribouTopology, HexahedronQuadraticFromIndices) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using namespace sofa::helper;
    using namespace sofa::core::objectmodel;

    using PointID  = sofa::core::topology::Topology::PointID;

    using Mesh = io::VTKReader<_3D, PointID>::MeshType;
    auto reader = io::VTKReader<_3D, PointID>::Read(executable_directory_path + "/meshes/3D_hexahedron_quadratic.vtk");
    using Domain = Mesh::Domain<Hexahedron20, PointID>;

    auto mesh = reader.mesh();
    EXPECT_EQ(mesh.number_of_nodes(), 425);
    EXPECT_EQ(mesh.number_of_domains(), 2);

    // First domain is the quad surface, second is the hexa volumetric domain
    // Get the second one
    const auto * domain = dynamic_cast<const Domain * >(mesh.domain(1));

    MessageDispatcher::addHandler( MainGtestMessageHandler::getInstance() ) ;
    EXPECT_MSG_NOEMIT(Error);

    setSimulation(new sofa::simulation::graph::DAGSimulation());
    auto root = getSimulation()->createNewNode("root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");

    // Add a mechanical object required for the creation of an internal mesh
    auto mo = dynamic_cast<sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3Types> *>(
            createObject(root, "MechanicalObject", {{"name", "mo"}}).get()
    );
    mo->resize(mesh.number_of_nodes());
    auto positions = mo->writePositions();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto p = mesh.position(i);
        positions[i] = {p[0], p[1], p[2]};
    }

    // Add the CaribouTopology component
    auto topo = dynamic_cast<SofaCaribou::topology::CaribouTopology<Hexahedron20> *> (
            createObject(root, "CaribouTopology", {{"template", "Hexahedron20"}}).get()
    );
    EXPECT_NE(topo, nullptr);

    // Set the indices
    using DataIndices = Data<sofa::type::vector<sofa::type::fixed_array<PointID, 20>>>;
    auto indices = WriteOnlyAccessor<DataIndices> (dynamic_cast<DataIndices*>(topo->findData("indices")));
    indices.resize(domain->number_of_elements());

    for (std::size_t element_id = 0; element_id < indices.size(); ++element_id) {
        auto & element_indices = indices[element_id];
        for (std::size_t node_id = 0; node_id < 20; ++node_id) {
            element_indices[node_id] = domain->element_indices(element_id)[node_id];
        }
    }

    getSimulation()->init(root.get());

    // Make sure the created internal domain correspond to the one in the mesh file
    const auto * hexa_domain = topo->domain();
    ASSERT_NE(hexa_domain, nullptr);

    FLOATING_POINT_TYPE volume = 0;
    for (UNSIGNED_INTEGER_TYPE hexa_id = 0; hexa_id < hexa_domain->number_of_elements(); ++hexa_id) {
        auto hexa = hexa_domain->element(hexa_id);
        for (const auto & g : hexa.gauss_nodes()) {
            const auto x = g.position;
            const auto w = g.weight;
            const auto J = hexa.jacobian(x);
            const auto detJ = J.determinant();
            volume += w*abs(detJ);
        }
    }
    double exact_volume = 10*10*10;
    EXPECT_LE( abs((volume-exact_volume)/exact_volume), 0.001); // 1% error max
}
