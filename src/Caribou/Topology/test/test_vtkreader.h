#ifndef CARIBOU_TOPOLOGY_TEST_VTKREADER_H
#define CARIBOU_TOPOLOGY_TEST_VTKREADER_H

#include <Caribou/Topology/UnstructuredMesh.h>
#include <Caribou/Topology/IO/VTKReader.h>

using namespace caribou::topology;

TEST(VTKReader, Read) {
    auto reader = io::VTKReader<3>::Read (executable_directory_path+"/meshes/1D_linear.vtk");
    auto mesh_1D_linear = reader.mesh();
    EXPECT_EQ(mesh_1D_linear.number_of_nodes(), 11);
    EXPECT_EQ(mesh_1D_linear.number_of_domains(), 1);
    EXPECT_EQ(mesh_1D_linear.domains()[0].first, "domain_1");
    EXPECT_EQ(mesh_1D_linear.domains()[0].second->number_of_elements(), 10);

}

#endif //CARIBOU_TOPOLOGY_TEST_VTKREADER_H
