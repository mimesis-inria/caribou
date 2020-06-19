#include <gtest/gtest.h>
#include "topology_test.h"
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/Domain.h>
#include <Caribou/Geometry/Segment.h>

TEST(Domain, Segment) {
    using namespace caribou;
    using namespace caribou::geometry;
    using namespace caribou::topology;

    using Mesh = Mesh<_3D>;

    { // Segment Linear (copy)
        using Segment = Segment<_3D, Linear>;
        using Domain = Mesh::Domain<Segment>;

        Mesh mesh;
        Domain::ElementsIndices indices(4, 2); // Four elements of two nodes each
        indices << 1, 2,
                   2, 3,
                   3, 4,
                   5, 6;

        Domain * domain = mesh.add_domain<Segment>(indices);
        EXPECT_MATRIX_EQUAL(domain->element_indices(0), indices.row(0).transpose());
        EXPECT_MATRIX_EQUAL(domain->element_indices(3), indices.row(3).transpose());
    }

    { // Segment Linear (reference)
        using Segment = Segment<_3D, Linear>;
        using Domain = Mesh::Domain<Segment>;

        Domain::ElementsIndices indices(4, 2); // Four elements of two nodes each
        indices << 1, 2,
            2, 3,
            3, 4,
            5, 6;

        Mesh mesh;
        Domain * domain = mesh.add_domain<Segment>(&indices);
        indices(0,0) = 0;
        indices(3,0) = 0;
        EXPECT_MATRIX_EQUAL(domain->element_indices(0), indices.row(0).transpose());
        EXPECT_MATRIX_EQUAL(domain->element_indices(3), indices.row(3).transpose());
    }

    { // Segment Linear (map)
        using Segment = Segment<_3D, Linear>;
        using Domain = Mesh::Domain<Segment>;

        UNSIGNED_INTEGER_TYPE indices[8] = {1, 2, 2, 3, 3, 4, 5, 6};
        Mesh mesh;
        Domain * domain = mesh.add_domain<Segment>(indices, 4, 2);
        indices[0] = 0;
        indices[7] = 0;
        EXPECT_EQ(domain->element_indices(0)[0], indices[0]);
        EXPECT_EQ(domain->element_indices(0)[1], indices[1]);
        EXPECT_EQ(domain->element_indices(3)[0], indices[6]);
        EXPECT_EQ(domain->element_indices(3)[1], indices[7]);
    }

    { // Segment Linear (map with stride)
        using Segment = Segment<_3D, Linear>;
        using Domain = Mesh::Domain<Segment>;

        // Two unused between columns, three unused betwen rows
        UNSIGNED_INTEGER_TYPE indices[25] = {1, 0, 0, 2, 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 3, 0, 0, 4, 0, 0, 0, 5, 0, 0, 6};
        Mesh mesh;
        Domain * domain = mesh.add_domain<Segment>(indices, 4, 2, 7, 3);
        indices[0] = 0;
        indices[24] = 0;
        EXPECT_EQ(domain->element_indices(0)[0], indices[0]);
        EXPECT_EQ(domain->element_indices(0)[1], indices[3]);
        EXPECT_EQ(domain->element_indices(3)[0], indices[21]);
        EXPECT_EQ(domain->element_indices(3)[1], indices[24]);
    }
}