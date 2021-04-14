project(Mechanics)

set(HEADER_FILES
    Elasticity/Strain.h)

#Eigen3
if(NOT ${Eigen3_FOUND})
    find_package(Eigen3 QUIET REQUIRED NO_MODULE)
endif()

add_library(${PROJECT_NAME} INTERFACE)
add_library(Caribou::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

# Add the target to the component of the same name. This will enable to do:
# CMakeLists.txt
#    find_package(Caribou COMPONENTS Mechanics)
#    add_library(my_cool_library SHARED)
#    target_link_libraries(my_cool_library PRIVATE Caribou::Mechanics)
caribou_add_target_to_component (
    TARGET_NAME    ${PROJECT_NAME}
    COMPONENT_NAME ${PROJECT_NAME}
    PUBLIC_HEADERS ${HEADER_FILES}
    PACKAGE_NAME   Caribou
)

# This is required to be able to do "find_package(Caribou COMPONENTS Mechanics)"
caribou_add_component_to_package (
    COMPONENT_NAME ${PROJECT_NAME}
    PACKAGE_NAME   Caribou
)
