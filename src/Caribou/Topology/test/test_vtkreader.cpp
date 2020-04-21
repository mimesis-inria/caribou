#include <gtest/gtest.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/IO/VTKReader.h>
#include "topology_test.h"

TEST(VTKReader, Segment) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;

    { // Linear
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/1D_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 11);
        EXPECT_EQ(mesh.number_of_domains(), 1);
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 10);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto domain = dynamic_cast<const Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
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
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/1D_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 21);
        EXPECT_EQ(mesh.number_of_domains(), 1);
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 10);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto domain = dynamic_cast<const Domain<Segment<_3D, Quadratic>> * >(mesh.domain(0));
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
        auto reader = io::VTKReader<_2D>::Read(executable_directory_path + "/meshes/2D_triangle_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 25);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto segment_domain = dynamic_cast<const Domain<Segment<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a segment to a triangle domain will return nullptr
        auto bad_domain = dynamic_cast<const Domain<Triangle<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain, nullptr);

        // Triangle domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 32);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        auto triangle_domain = dynamic_cast<const Domain<Triangle<_2D, Linear>> * >(mesh.domain(1));
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
        auto reader = io::VTKReader<_2D>::Read(executable_directory_path + "/meshes/2D_triangle_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 81);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto segment_domain = dynamic_cast<const Domain<Segment<_2D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a quadratic segment to a linear segment domain will return nullptr
        auto bad_domain1 = dynamic_cast<const Domain<Segment<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a segment to a triangle domain will return nullptr
        auto bad_domain2 = dynamic_cast<const Domain<Triangle<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain2, nullptr);

        // Triangle domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 32);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        auto triangle_domain = dynamic_cast<const Domain<Triangle<_2D, Quadratic>> * >(mesh.domain(1));
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
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_triangle_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 25);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto segment_domain = dynamic_cast<const Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a segment to a triangle domain will return nullptr
        auto bad_domain1 = dynamic_cast<const Domain<Triangle<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a 3D domain to a 2D domain will return nullptr
        auto bad_domain2 = dynamic_cast<const Domain<Triangle<_2D, Linear>> * >(mesh.domain(1));
        ASSERT_EQ(bad_domain2, nullptr);

        // Triangle domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 32);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        auto triangle_domain = dynamic_cast<const Domain<Triangle<_3D, Linear>> * >(mesh.domain(1));
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
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/2D_triangle_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 81);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto segment_domain = dynamic_cast<const Domain<Segment<_3D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a quadratic segment to a linear segment domain will return nullptr
        auto bad_domain1 = dynamic_cast<const Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a segment to a triangle domain will return nullptr
        auto bad_domain2 = dynamic_cast<const Domain<Triangle<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain2, nullptr);

        // Triangle domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 32);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        auto triangle_domain = dynamic_cast<const Domain<Triangle<_3D, Quadratic>> * >(mesh.domain(1));
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
        auto reader = io::VTKReader<_2D>::Read(executable_directory_path + "/meshes/2D_quad_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 25);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto segment_domain = dynamic_cast<const Domain<Segment<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a segment to a quad domain will return nullptr
        auto bad_domain = dynamic_cast<const Domain<Quad<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain, nullptr);

        // Quad domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        auto quad_domain = dynamic_cast<const Domain<Quad<_2D, Linear>> * >(mesh.domain(1));
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
        auto reader = io::VTKReader<_2D>::Read(executable_directory_path + "/meshes/2D_quad_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 65);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto segment_domain = dynamic_cast<const Domain<Segment<_2D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);


        // Trying to cast a quadratic segment to a linear segment domain will return nullptr
        auto bad_domain1 = dynamic_cast<const Domain<Segment<_2D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a quad to a triangle domain will return nullptr
        auto bad_domain2 = dynamic_cast<const Domain<Triangle<_2D, Quadratic>> * >(mesh.domain(1));
        ASSERT_EQ(bad_domain2, nullptr);

        // Quad domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        auto quad_domain = dynamic_cast<const Domain<Quad<_2D, Quadratic>> * >(mesh.domain(1));
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
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_quad_linear.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 25);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto segment_domain = dynamic_cast<const Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);

        // Trying to cast a segment to a quad domain will return nullptr
        auto bad_domain = dynamic_cast<const Domain<Quad<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain, nullptr);

        // Quad domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        auto quad_domain = dynamic_cast<const Domain<Quad<_3D, Linear>> * >(mesh.domain(1));
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
        auto reader = io::VTKReader<_3D>::Read(executable_directory_path + "/meshes/3D_quad_quadratic.vtk");
        auto mesh = reader.mesh();
        EXPECT_EQ(mesh.number_of_nodes(), 65);
        EXPECT_EQ(mesh.number_of_domains(), 2);

        // Contour segments
        EXPECT_EQ(mesh.domains()[0].first, "domain_1");
        EXPECT_EQ(mesh.domain(0)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(0), mesh.domain("domain_1"));

        auto segment_domain = dynamic_cast<const Domain<Segment<_3D, Quadratic>> * >(mesh.domain(0));
        ASSERT_NE(segment_domain, nullptr);


        // Trying to cast a quadratic segment to a linear segment domain will return nullptr
        auto bad_domain1 = dynamic_cast<const Domain<Segment<_3D, Linear>> * >(mesh.domain(0));
        ASSERT_EQ(bad_domain1, nullptr);

        // Trying to cast a quad to a triangle domain will return nullptr
        auto bad_domain2 = dynamic_cast<const Domain<Triangle<_3D, Quadratic>> * >(mesh.domain(1));
        ASSERT_EQ(bad_domain2, nullptr);

        // Quad domain
        EXPECT_EQ(mesh.domains()[1].first, "domain_2");
        EXPECT_EQ(mesh.domain(1)->number_of_elements(), 16);
        EXPECT_EQ(mesh.domain(1), mesh.domain("domain_2"));
        auto quad_domain = dynamic_cast<const Domain<Quad<_3D, Quadratic>> * >(mesh.domain(1));
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