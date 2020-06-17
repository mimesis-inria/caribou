#include <gtest/gtest.h>
#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/IO/VTKReader.h>
#include "topology_test.h"

TEST(VTKReader, Segment) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

    { // Linear
        using Mesh = io::VTKReader<_3D>::MeshType;
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/1D_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 11);
        EXPECT_EQ(mesh.number_of_domains(), 1);
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 10);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * domain = dynamic_cast<const Mesh::Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(domain, nullptr);

        FLOATING_POINT_TYPE segment_length = 0;
        for (UNSIGNED_INTEGER_TYPE segment_id = 0; segment_id < domain->number_of_elements(); ++segment_id) {
            auto segment = domain->element(segment_id);
            for (const auto & g : segment.gauss_nodes()) {
                const auto x = g.position;
                const auto w = g.weight;
                const auto J =  segment.jacobian(x);
                const auto detJ = J.norm();
                segment_length += w*detJ;
            }
        }
        EXPECT_FLOAT_EQ(segment_length, 10.);
    }

    { // Quadratic
        using Mesh = io::VTKReader<_3D>::MeshType;
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/1D_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 21);
        EXPECT_EQ(mesh.number_of_domains(), 1);
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 10);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * domain = dynamic_cast<const Mesh::Domain<Segment<_3D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(domain, nullptr);

        FLOATING_POINT_TYPE segment_length = 0;
        for (UNSIGNED_INTEGER_TYPE segment_id = 0; segment_id < domain->number_of_elements(); ++segment_id) {
            auto segment = domain->element(segment_id);
            for (const auto & g : segment.gauss_nodes()) {
                const auto x = g.position;
                const auto w = g.weight;
                const auto J =  segment.jacobian(x);
                const auto detJ = J.norm();
                segment_length += w*detJ;
            }
        }
        EXPECT_FLOAT_EQ(segment_length, 10.);
    }
}

TEST(VTKReader, Triangle2D) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

    { // Linear 2D
        using Mesh = io::VTKReader<_2D>::MeshType;
        auto reader = io::VTKReader<_2D>::Read(executable_directory_path + "/meshes/2D_triangle_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 25);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * segment_domain = dynamic_cast<const Mesh::Domain<Segment<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a segment to a triangle domain will return nullptr
        const auto * bad_domain = dynamic_cast<const Mesh::Domain<Triangle<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain, nullptr);

        // Triangle domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 32);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * triangle_domain = dynamic_cast<const Mesh::Domain<Triangle<_2D, Linear>> * >(mesh.domain(1));
        ASSERT_NE(triangle_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), triangle_domain->number_of_elements());

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

    { // Quadratic 2D
        using Mesh = io::VTKReader<_2D>::MeshType;
        auto reader = io::VTKReader<_2D>::Read(executable_directory_path + "/meshes/2D_triangle_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 81);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * segment_domain = dynamic_cast<const Mesh::Domain<Segment<_2D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a quadratic segment to a linear segment domain will return nullptr
        const auto *  bad_domain1 = dynamic_cast<const Mesh::Domain<Segment<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a segment to a triangle domain will return nullptr
        const auto *  bad_domain2 = dynamic_cast<const Mesh::Domain<Triangle<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain2, nullptr);

        // Triangle domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 32);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto *  triangle_domain = dynamic_cast<const Mesh::Domain<Triangle<_2D, Quadratic>> * >(mesh.domain(1));
        ASSERT_NE(triangle_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), triangle_domain->number_of_elements());

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
}

TEST(VTKReader, Triangle3D) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

    { // Linear 3D
        using Mesh = io::VTKReader<_3D>::MeshType;
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_triangle_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 25);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * segment_domain = dynamic_cast<const Mesh::Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a segment to a triangle domain will return nullptr
        const auto * bad_domain1 = dynamic_cast<const Mesh::Domain<Triangle<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a 3D domain to a 2D domain will return nullptr
        const auto * bad_domain2 = dynamic_cast<const Mesh::Domain<Triangle<_2D, Linear>> * >(mesh.domain(1));
        ASSERT_EQ(bad_domain2, nullptr);

        // Triangle domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 32);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * triangle_domain = dynamic_cast<const Mesh::Domain<Triangle<_3D, Linear>> * >(mesh.domain(1));
        ASSERT_NE(triangle_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), triangle_domain->number_of_elements());

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

    { // Quadratic 3D
        using Mesh = io::VTKReader<_3D>::MeshType;
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/2D_triangle_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 81);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * segment_domain = dynamic_cast<const Mesh::Domain<Segment<_3D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a quadratic segment to a linear segment domain will return nullptr
        const auto * bad_domain1 = dynamic_cast<const Mesh::Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a segment to a triangle domain will return nullptr
        const auto * bad_domain2 = dynamic_cast<const Mesh::Domain<Triangle<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain2, nullptr);

        // Triangle domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 32);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * triangle_domain = dynamic_cast<const Mesh::Domain<Triangle<_3D, Quadratic>> * >(mesh.domain(1));
        ASSERT_NE(triangle_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), triangle_domain->number_of_elements());

        FLOATING_POINT_TYPE area = 0;
        for (UNSIGNED_INTEGER_TYPE triangle_id = 0; triangle_id < triangle_domain->number_of_elements(); ++triangle_id) {
            auto triangle = triangle_domain->element(triangle_id);
            for (const auto & g : triangle.gauss_nodes()) {
                const auto x = g.position;
                const auto w = g.weight;
                const auto J = triangle.jacobian(x);
                const auto detJ = J.col(0).cross(J.col(1)).norm();
                area += w*abs(detJ);
            }
        }
        EXPECT_NEAR(area, 100., 1e-3);
    }
}

TEST(VTKReader, Quad2D) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

    { // Linear 2D
        using Mesh = io::VTKReader<_2D>::MeshType;
        auto reader = io::VTKReader<_2D>::Read(executable_directory_path + "/meshes/2D_quad_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 25);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * segment_domain = dynamic_cast<const Mesh::Domain<Segment<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a segment to a quad domain will return nullptr
        const auto * bad_domain = dynamic_cast<const Mesh::Domain<Quad<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain, nullptr);

        // Quad domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * quad_domain = dynamic_cast<const Mesh::Domain<Quad<_2D, Linear>> * >(mesh.domain(1));
        ASSERT_NE(quad_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), quad_domain->number_of_elements());

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

    { // Quadratic 2D
        using Mesh = io::VTKReader<_2D>::MeshType;
        auto reader = io::VTKReader<_2D>::Read(executable_directory_path + "/meshes/2D_quad_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 65);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * segment_domain = dynamic_cast<const Mesh::Domain<Segment<_2D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);


        // Trying to cast a quadratic segment to a linear segment domain will return nullptr
        const auto * bad_domain1 = dynamic_cast<const Mesh::Domain<Segment<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a quad to a triangle domain will return nullptr
        const auto * bad_domain2 = dynamic_cast<const Mesh::Domain<Triangle<_2D, Quadratic>> * >(mesh.domain(1));
        ASSERT_EQ(bad_domain2, nullptr);

        // Quad domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * quad_domain = dynamic_cast<const Mesh::Domain<Quad<_2D, Quadratic>> * >(mesh.domain(1));
        ASSERT_NE(quad_domain, nullptr);
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), quad_domain->number_of_elements());

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
}

TEST(VTKReader, Quad3D) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

    { // Linear 3D
        using Mesh = io::VTKReader<_3D>::MeshType;
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_quad_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 25);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * segment_domain = dynamic_cast<const Mesh::Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a segment to a quad domain will return nullptr
        const auto * bad_domain = dynamic_cast<const Mesh::Domain<Quad<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain, nullptr);

        // Quad domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * quad_domain = dynamic_cast<const Mesh::Domain<Quad<_3D, Linear>> * >(mesh.domain(1));
        ASSERT_NE(quad_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), quad_domain->number_of_elements());

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

    { // Quadratic 3D
        using Mesh = io::VTKReader<_3D>::MeshType;
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_quad_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 65);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * segment_domain = dynamic_cast<const Mesh::Domain<Segment<_3D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);


        // Trying to cast a quadratic segment to a linear segment domain will return nullptr
        const auto * bad_domain1 = dynamic_cast<const Mesh::Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a quad to a triangle domain will return nullptr
        const auto * bad_domain2 = dynamic_cast<const Mesh::Domain<Triangle<_3D, Quadratic>> * >(mesh.domain(1));
        ASSERT_EQ(bad_domain2, nullptr);

        // Quad domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * quad_domain = dynamic_cast<const Mesh::Domain<Quad<_3D, Quadratic>> * >(mesh.domain(1));
        ASSERT_NE(quad_domain, nullptr);
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), quad_domain->number_of_elements());

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
}

TEST(VTKReader, Tetrahedron) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

    double pi = 3.14159265358979323846;
    double r = 5;
    using Mesh = io::VTKReader<_3D>::MeshType;

    { // Linear
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_tetrahedron_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 206);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Surface triangles
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 316);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * triangle_domain = dynamic_cast<const Mesh::Domain<Triangle<_3D, Linear>> * >(mesh.domain(0));
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
        double exact_area = 4*pi*r*r;
        EXPECT_LE( abs((area-exact_area)/exact_area), 0.1);

        // Trying to cast a triangle to a tetra domain will return nullptr
        const auto * bad_domain = dynamic_cast<const Mesh::Domain<Tetrahedron<Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain, nullptr);

        // Tetra domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 681);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * tetra_domain = dynamic_cast<const Mesh::Domain<Tetrahedron<Linear>> * >(mesh.domain(1));
        ASSERT_NE(tetra_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), tetra_domain->number_of_elements());

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
        EXPECT_LE( abs((volume-exact_volume)/exact_volume), 0.1);
    }

    { // Quadratic
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_tetrahedron_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 1250);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Surface triangles
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 316);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * triangle_domain = dynamic_cast<const Mesh::Domain<Triangle<_3D, Quadratic>> * >(mesh.domain(0));
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
        double exact_area = 4*pi*r*r;
        EXPECT_LE( abs((area-exact_area)/exact_area), 0.0005);

        // Trying to cast a triangle to a tetra domain will return nullptr
        const auto * bad_domain1 = dynamic_cast<const Mesh::Domain<Tetrahedron<Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a quadratic triangle to a linear triangle domain will return nullptr
        const auto * bad_domain2 = dynamic_cast<const Mesh::Domain<Triangle<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain2, nullptr);

        // Tetra domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 681);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * tetra_domain = dynamic_cast<const Mesh::Domain<Tetrahedron<Quadratic>> * >(mesh.domain(1));
        ASSERT_NE(tetra_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), tetra_domain->number_of_elements());

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
        EXPECT_LE( abs((volume-exact_volume)/exact_volume), 0.0005);
    }
}

TEST(VTKReader, Hexahedron) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

    using Mesh = io::VTKReader<_3D>::MeshType;

    { // Linear
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_hexahedron_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 125);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // 6 surface quads
        FLOATING_POINT_TYPE area = 0;
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16*6);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * quad_domain = dynamic_cast<const Mesh::Domain<Quad<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(quad_domain, nullptr);

        for (UNSIGNED_INTEGER_TYPE quad_id = 0; quad_id < quad_domain->number_of_elements(); ++quad_id) {
            auto quad = quad_domain->element(quad_id);
            for (const auto & g : quad.gauss_nodes()) {
                const auto x = g.position;
                const auto w = g.weight;
                const auto J = quad.jacobian(x);
                const auto detJ = J.col(0).cross(J.col(1)).norm();
                area += w*detJ;
            }
        }
        double exact_area = 6*100;
        EXPECT_LE( abs((area-exact_area)/exact_area), 0.001);

        // Trying to cast an hexa to a tetra domain will return nullptr
        const auto * bad_domain = dynamic_cast<const Mesh::Domain<Tetrahedron<Linear>> * >(mesh.domain(1));
        ASSERT_EQ(bad_domain, nullptr);

        // Hexa domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 64);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * hexa_domain = dynamic_cast<const Mesh::Domain<Hexahedron<Linear>> * >(mesh.domain(1));
        ASSERT_NE(hexa_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), hexa_domain->number_of_elements());

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
        EXPECT_LE( abs((volume-exact_volume)/exact_volume), 0.001);
    }

    { // Quadratic
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_hexahedron_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 425);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // 6 surface quads
        FLOATING_POINT_TYPE area = 0;
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16*6);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        const auto * quad_domain = dynamic_cast<const Mesh::Domain<Quad<_3D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(quad_domain, nullptr);

        for (UNSIGNED_INTEGER_TYPE quad_id = 0; quad_id < quad_domain->number_of_elements(); ++quad_id) {
            auto quad = quad_domain->element(quad_id);
            for (const auto & g : quad.gauss_nodes()) {
                const auto x = g.position;
                const auto w = g.weight;
                const auto J = quad.jacobian(x);
                const auto detJ = J.col(0).cross(J.col(1)).norm();
                area += w*detJ;
            }
        }
        double exact_area = 6*100;
        EXPECT_LE( abs((area-exact_area)/exact_area), 0.001);

        // Trying to cast an hexa to a tetra domain will return nullptr
        const auto * bad_domain = dynamic_cast<const Mesh::Domain<Tetrahedron<Linear>> * >(mesh.domain(1));
        ASSERT_EQ(bad_domain, nullptr);

        // Hexa domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 64);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        const auto * hexa_domain = dynamic_cast<const Mesh::Domain<Hexahedron<Quadratic>> * >(mesh.domain(1));
        ASSERT_NE(hexa_domain, nullptr);
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), hexa_domain->number_of_elements());

        FLOATING_POINT_TYPE volume = 0;
        for (UNSIGNED_INTEGER_TYPE hexa_id = 0; hexa_id < hexa_domain->number_of_elements(); ++hexa_id) {
            auto hexa = hexa_domain->element(hexa_id);
            for (const auto & g : hexa.gauss_nodes()) {
                const auto x = g.position;
                const auto w = g.weight;
                const auto J = hexa.jacobian(x);
                const auto detJ = J.determinant();
                volume += w*detJ;
            }
        }
        double exact_volume = 10*10*10;
        EXPECT_LE( abs((volume-exact_volume)/exact_volume), 0.001);
    }
}