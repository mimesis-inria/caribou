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
#include <sofa/version.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION >= 220600)
#include <sofa/simulation/graph/init.h>
#else
#include <SofaSimulationGraph/init.h>
#endif
#if (defined(SOFA_VERSION) && SOFA_VERSION >= 201200)
#include <SofaBaseMechanics/initSofaBaseMechanics.h>
#include <SofaBaseUtils/initSofaBaseUtils.h>
#else
#include <SofaBaseMechanics/initBaseMechanics.h>
#include <SofaBaseUtils/initBaseUtils.h>
#endif
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
    sofa::component::initBaseMechanics();
    sofa::component::initBaseUtils();
#endif

    int ret = RUN_ALL_TESTS();
    sofa::simulation::graph::cleanup();
    return ret;
}
