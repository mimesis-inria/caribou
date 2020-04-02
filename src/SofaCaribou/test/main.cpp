#include <gtest/gtest.h>
#include <SofaSimulationGraph/init.h>

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    sofa::simulation::graph::init();
    int ret = RUN_ALL_TESTS();
    sofa::simulation::graph::cleanup();
    return ret;
}
