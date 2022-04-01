#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/FictitiousGrid.h>

DISABLE_ALL_WARNINGS_BEGIN
#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
#include <sofa/helper/testing/BaseTest.h>
#else
#include <sofa/testing/BaseTest.h>
#endif
#include <sofa/simulation/Node.h>
#include <SofaSimulationGraph/DAGSimulation.h>
#include <SofaSimulationGraph/SimpleApi.h>
#include <sofa/helper/system/PluginManager.h>
DISABLE_ALL_WARNINGS_END

#include "../sofacaribou_test.h"

using sofa::helper::system::PluginManager ;
using namespace sofa::simulation;
using namespace sofa::simpleapi;

#if (defined(SOFA_VERSION) && SOFA_VERSION < 210600)
using namespace sofa::helper::testing;
#else
using namespace sofa::testing;
#endif

class FictitiousGrid : public BaseTest {
    void SetUp() override {
        setSimulation(new sofa::simulation::graph::DAGSimulation()) ;
        root = getSimulation()->createNewNode("root");
#if (defined(SOFA_VERSION) && SOFA_VERSION >= 201200)
        createObject(root, "RequiredPlugin", {{"pluginName", "SofaGeneralLoader"}});
#endif
        createObject(root, "DefaultAnimationLoop");
        createObject(root, "DefaultVisualManagerLoop");
    }
    void TearDown() override {
        root.reset();
        setSimulation(nullptr);
    }

protected:
    sofa::simulation::Node::SPtr root;
};

TEST_F(FictitiousGrid, Liver) {
    EXPECT_MSG_NOEMIT(Error, Warning) ;
    createObject(root, "MeshSTLLoader", {{"name", "loader"}, {"filename", executable_directory_path + "/meshes/deformed_liver_surface.stl"}});
    auto grid = dynamic_cast<SofaCaribou::topology::FictitiousGrid<sofa::defaulttype::Vec3Types> *>(createObject(root, "FictitiousGrid", {
        {"printLog", "0"},
        {"n", "15 15 15"},
        {"maximum_number_of_subdivision_levels", std::to_string(4)},
        {"surface_positions", "@./loader.position"},
        {"surface_triangles", "@./loader.triangles"}
    }).get());

    getSimulation()->init(root.get());

    // Get the volume
    FLOATING_POINT_TYPE volume = 0.;
    for (std::size_t element_id = 0; element_id < grid->number_of_cells(); ++element_id) {
        const auto gauss_nodes = grid->get_gauss_nodes_of_cell(element_id, 0);
        for (const auto & gauss_node : gauss_nodes) {
            volume += gauss_node.second;
        }
    }

    EXPECT_NEAR(volume, 3414171, 5);
}