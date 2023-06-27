#ifdef LEGACY_CXX
#include <experimental/filesystem>
namespace fs = ::std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = ::std::filesystem;
#endif

#include <gtest/gtest.h>

#include <SofaCaribou/config.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/simulation/graph/init.h>
#include <sofa/component/mass/init.h>
#include <sofa/component/statecontainer/init.h>
DISABLE_ALL_WARNINGS_END

#include "sofacaribou_test.h"

std::string executable_directory_path;

int main(int argc, char **argv) {
#ifdef LEGACY_CXX
    executable_directory_path = fs::canonical(fs::path(argv[0])).parent_path();
#else
    executable_directory_path = weakly_canonical(fs::path(argv[0])).parent_path().string();
#endif
    testing::InitGoogleTest(&argc, argv);
    sofa::simulation::graph::init();
#if (defined(SOFA_VERSION) && SOFA_VERSION >= 201200)
    sofa::component::initSofaBaseMechanics();
    sofa::component::initSofaBaseUtils();
#else
    sofa::component::mass::init();
    sofa::component::statecontainer::init();
#endif

    int ret = RUN_ALL_TESTS();
    sofa::simulation::graph::cleanup();
    return ret;
}
