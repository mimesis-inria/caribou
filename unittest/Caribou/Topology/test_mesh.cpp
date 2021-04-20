#ifndef CARIBOU_TOPOLOGY_TEST_MESH_H
#define CARIBOU_TOPOLOGY_TEST_MESH_H

#include <gtest/gtest.h>
#include "topology_test.h"
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/Domain.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Tetrahedron.h>

using namespace caribou;
using namespace caribou::geometry;
using caribou::topology::Mesh;

TEST(Mesh, Constructor_std_copy) {
    // Construct by copy (from std::vector)
    std::vector<Mesh<3>::WorldCoordinates> initial_positions = {{0,0,0}, {1,1,1}};
    Mesh<3> mesh (initial_positions);
    auto positions = mesh.positions({1, 0});
    EXPECT_MATRIX_EQUAL(mesh.position(1), positions.row(0));
    EXPECT_MATRIX_EQUAL(mesh.position(0), positions.row(1));
    EXPECT_NE(initial_positions.data()->data(), mesh.position(0).data());
    EXPECT_NE(initial_positions[1].data(), mesh.position(1).data());
    EXPECT_NE(&initial_positions[1][2], &mesh.position(1).coeff(2));
}

TEST(Mesh, Constructor_eigen_copy) {
    // Construct by copy (from Eigen::Matrix)
    Eigen::Matrix<FLOATING_POINT_TYPE, 2, 3> initial_positions;
    initial_positions << 0, 0, 0,
                         1, 1, 1;
    Mesh mesh (initial_positions);
    auto positions = mesh.positions({1, 0});
    EXPECT_MATRIX_EQUAL(mesh.position(1), positions.row(0));
    EXPECT_MATRIX_EQUAL(mesh.position(0), positions.row(1));
    EXPECT_NE(initial_positions.data(), mesh.position(0).data());
    EXPECT_NE(&initial_positions(1,0), mesh.position(1).data());
    EXPECT_NE(&initial_positions(1,2), &mesh.position(1).coeff(2));

    // Construct by copy (from Eigen::Ref)
    Eigen::Ref<Eigen::Matrix<FLOATING_POINT_TYPE, 2, 3>> pos_ref(initial_positions);
    Mesh mesh_ref (pos_ref);
    positions = mesh_ref.positions({1, 0});
    EXPECT_MATRIX_EQUAL(mesh_ref.position(0), initial_positions.row(0));
    EXPECT_MATRIX_EQUAL(mesh_ref.position(1), initial_positions.row(1));
    EXPECT_MATRIX_EQUAL(mesh_ref.position(1), positions.row(0));
    EXPECT_MATRIX_EQUAL(mesh_ref.position(0), positions.row(1));
    EXPECT_NE(initial_positions.data(), mesh_ref.position(0).data());
    EXPECT_NE(&initial_positions(1,0), mesh_ref.position(1).data());
    EXPECT_NE(&initial_positions(1,2), &mesh_ref.position(1).coeff(2));

    // Copy constructor
    auto mesh2 = mesh;
    positions = mesh2.positions({1, 0});
    EXPECT_MATRIX_EQUAL(mesh2.position(0), initial_positions.row(0));
    EXPECT_MATRIX_EQUAL(mesh2.position(1), initial_positions.row(1));
    EXPECT_MATRIX_EQUAL(mesh2.position(1), positions.row(0));
    EXPECT_MATRIX_EQUAL(mesh2.position(0), positions.row(1));
    EXPECT_NE(initial_positions.data(), mesh2.position(0).data());
    EXPECT_NE(&initial_positions(1,0), mesh2.position(1).data());
    EXPECT_NE(&initial_positions(1,2), &mesh2.position(1).coeff(2));
}

TEST(Mesh, Segment) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using Mesh = Mesh<_3D>;
    using Domain = Mesh::Domain<Segment<_3D>>;

    Mesh mesh({{0, 0, 0}, {1, 1, 1}});

    Domain::ElementsIndices domain_indices(1, 2); // One element of two nodes
    domain_indices << 0, 1;

    Domain * domain = mesh.add_domain<Segment<3>>("segments", domain_indices);

    EXPECT_EQ(domain->number_of_elements(), 1);

    {
        const auto &indices = domain->element_indices(0);
        EXPECT_MATRIX_EQUAL(indices, Domain::ElementIndices ({0, 1}));
    }

    EXPECT_ANY_THROW(mesh.add_domain<Segment<3>>("segments", domain_indices));

    {
        auto positions = mesh.positions(domain->element_indices(0));
        const Segment<3> segment(positions);
        EXPECT_MATRIX_EQUAL(segment.center(), Segment<3>::WorldCoordinates(0.5, 0.5, 0.5));
    }

    mesh.remove(domain);
    EXPECT_EQ(mesh.number_of_domains(), 0);
}

TEST(Mesh, Triangle) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using Mesh = Mesh<_3D>;
    using Domain = Mesh::Domain<Triangle<_3D>>;

    Mesh mesh ({{50,50,33}, {60, 50, 21}, {55, 55, -4}});

    Domain::ElementsIndices domain_indices(1, 3); // One element of three nodes
    domain_indices << 2, 1, 0;

    Domain * domain = mesh.add_domain<Triangle<_3D, Linear>>("triangles", domain_indices);

    EXPECT_EQ(domain->number_of_elements(), 1);

    {
        const auto &indices = domain->element_indices(0);
        EXPECT_MATRIX_EQUAL(indices, Domain::ElementIndices({2, 1, 0}));
    }

    EXPECT_ANY_THROW(mesh.add_domain<Triangle<3>>("triangles", domain_indices));


    {
        auto positions = mesh.positions(domain->element_indices(0));
        const Triangle<3> triangle (positions);
        EXPECT_MATRIX_EQUAL(triangle.node(0), Triangle<_3D>::WorldCoordinates(55, 55, -4));
    }

    {
        auto triangle = domain->element(0);
        EXPECT_MATRIX_EQUAL(triangle.node(0), Triangle<_3D>::WorldCoordinates(55, 55, -4));
        EXPECT_MATRIX_EQUAL(triangle.node(1), Triangle<_3D>::WorldCoordinates(60, 50, 21));
        EXPECT_MATRIX_EQUAL(triangle.node(2), Triangle<_3D>::WorldCoordinates(50,50,33));
    }

    mesh.remove(domain);
    EXPECT_EQ(mesh.number_of_domains(), 0);
}

TEST(Mesh, Tetrahedron) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using Mesh = Mesh<_3D>;
    using Domain = Mesh::Domain<Tetrahedron<Linear>>;

    Mesh mesh ({{50,50,33}, {50,50,33}, {55, 55, -4}, {55, 55, -40}});

    Domain::ElementsIndices domain_indices(1, 4); // One element of four nodes
    domain_indices << 3, 2, 1, 0;

    Domain * domain = mesh.add_domain<Tetrahedron<Linear>>("tetras", domain_indices);

    EXPECT_EQ(domain->number_of_elements(), 1);

    {
        const auto &indices = domain->element_indices(0);
        EXPECT_MATRIX_EQUAL(indices, Domain::ElementIndices({3, 2, 1, 0}));
    }

    EXPECT_ANY_THROW(mesh.add_domain<Tetrahedron<Linear>>("tetras", domain_indices));


    {
        auto positions = mesh.positions(domain->element_indices(0));
        const Tetrahedron<Linear> tetra (positions);
        EXPECT_MATRIX_EQUAL(tetra.node(0), Tetrahedron<Linear>::WorldCoordinates(55, 55, -40));
    }

    {
        auto tetra = domain->element(0);
        EXPECT_MATRIX_EQUAL(tetra.node(1), Tetrahedron<Linear>::WorldCoordinates(55, 55, -4));
    }

    mesh.remove(domain);
    EXPECT_EQ(mesh.number_of_domains(), 0);
}

#endif //CARIBOU_TOPOLOGY_TEST_MESH_H