#include <filesystem>
#include <gtest/gtest.h>

#include "topology_test.h"

std::string executable_directory_path;

int main(int argc, char **argv) {
    executable_directory_path = weakly_canonical(std::filesystem::path(argv[0])).parent_path();
    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
