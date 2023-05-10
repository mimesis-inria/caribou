#include <SofaCaribou/config.h>
#include <SofaCaribou/Forcefield/TractionForcefield[Quad].h>
#include <SofaCaribou/Forcefield/TractionForcefield[Triangle].h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/testing/BaseTest.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/graph/DAGSimulation.h>
#include <sofa/simulation/graph/SimpleApi.h>
#include <sofa/helper/system/PluginManager.h>
DISABLE_ALL_WARNINGS_END

using sofa::helper::system::PluginManager ;
using namespace sofa::simulation;
using namespace sofa::simpleapi;
using namespace caribou;
using namespace caribou::geometry;
using namespace sofa::testing;

class TractionForcefield : public BaseTest {
    void SetUp() override {
        setSimulation(new sofa::simulation::graph::DAGSimulation()) ;
        root = getSimulation()->createNewNode("TractionForce");
    }
    void TearDown() override {
        root.reset();
        setSimulation(nullptr);
    }

protected:
    sofa::simulation::Node::SPtr root;
};

TEST_F(TractionForcefield, Triangle) {
    EXPECT_MSG_NOEMIT(Error, Warning) ;
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");
    createObject(root, "MechanicalObject", {{"position", "-1 0 1  1 0 1  -1 0 -1  1 0 -1  0 0 1  0 0 -1  -1 0 0  1 0 0  0 0 0"}});
    createObject(root, "TriangleSetTopologyContainer", {{"triangles", "7 5 8  8 2 6  4 6 0  1 8 4  7 3 5  8 5 2  4 8 6  1 7 8"}});
    auto traction = dynamic_cast<const SofaCaribou::forcefield::TractionForcefield<Triangle<_3D>> *>(
        createObject(root, "TractionForcefield", {{"traction", "0 5 0"}, {"slope", std::to_string(1/5.)}}).get()
    );

    EXPECT_NE(traction, nullptr);

    getSimulation()->init(root.get());
    auto total_load = dynamic_cast<sofa::core::objectmodel::Data<double> *>(traction->findData("total_load"));
    for (unsigned int step = 1; step <= 5; ++step) {
        getSimulation()->animate(root.get(), 1);
        EXPECT_DOUBLE_EQ(total_load->getValue(), 4*step) << "Total load at time step " << step << " is incorrect.";
    }
}

TEST_F(TractionForcefield, Quad) {
    EXPECT_MSG_NOEMIT(Error, Warning) ;
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");
    createObject(root, "MechanicalObject", {{"position", "-1 0 1  1 0 1  -1 0 -1  1 0 -1  0 0 1  0 0 -1  -1 0 0  1 0 0  0 0 0"}});
    createObject(root, "QuadSetTopologyContainer", {{"quads", "8 7 3 5  6 8 5 2  0 4 8 6  4 1 7 8"}});
    auto traction = createObject(root, "TractionForcefield", {{"traction", "0 5 0"}, {"slope", std::to_string(1/5.)}});

    getSimulation()->init(root.get());
    auto total_load = dynamic_cast<sofa::core::objectmodel::Data<double> *>(traction->findData("total_load"));
    for (unsigned int step = 1; step <= 5; ++step) {
        getSimulation()->animate(root.get(), 1);
        EXPECT_DOUBLE_EQ(total_load->getValue(), 4*step) << "Total load at time step " << step << " is incorrect.";
    }
}