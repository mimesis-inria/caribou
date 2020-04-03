#include <gtest/gtest.h>
#include <Caribou/Geometry/Segment.h>
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