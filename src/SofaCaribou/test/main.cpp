#ifdef LEGACY_CXX
#include <experimental/filesystem>
namespace fs = ::std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = ::std::filesystem;
#endif

#include <gtest/gtest.h>
#include <SofaSimulationGraph/init.h>
#include <SofaBaseMechanics/initSofaBaseMechanics.h>

#include "sofacaribou_test.h"

std::string executable_directory_path;

int main(int argc, char **argv) {
#ifdef LEGACY_CXX
    executable_directory_path = fs::canonical(fs::path(argv[0])).parent_path();
#else
    executable_directory_path = weakly_canonical(fs::path(argv[0])).parent_path();
#endif
    testing::InitGoogleTest(&argc, argv);
    sofa::simulation::graph::init();
    sofa::component::initSofaBaseMechanics();
    int ret = RUN_ALL_TESTS();
    sofa::simulation::graph::cleanup();
    return ret;
}
