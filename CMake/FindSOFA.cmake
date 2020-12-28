################################################################################
# Find SOFA framework libraries
#
#  SOFA_FOUND - Found SOFA libraries
#  SOFA_MODULE_DIR - Path to SOFA CMake modules
#  SOFA_ROOT - Path to SOFA installation
#
#  The environment variable SOFA_ROOT can be used to find SOFA.
#
#  Example usage:
#
#  find_package(SOFA COMPONENTS SofaCommon REQUIRED)
#  if(SOFA_FOUND)
#    target_link_libraries(${PROJECT_NAME} SofaEngine)
#  endif()
################################################################################

# Try to find SOFA's cmake module directory
find_package(SofaFramework CONFIG QUIET) # This defines SOFA_ROOT if SofaFrameworkConfig.cmake is found
find_path (
    SOFA_MODULE_DIR SofaFramework/SofaFrameworkConfig.cmake
    PATHS
        ${SOFA_ROOT}/lib/cmake
        ${SOFA_ROOT}/install/lib/cmake
        ${SOFA_ROOT}/build/install/lib/cmake
        ${SOFA_ROOT}/build/master/install/lib/cmake
        $ENV{SOFA_ROOT}/lib/cmake
        $ENV{SOFA_ROOT}/install/lib/cmake
        $ENV{SOFA_ROOT}/build/install/lib/cmake
        $ENV{SOFA_ROOT}/build/master/install/lib/cmake
)

if (SOFA_MODULE_DIR)
    get_filename_component(SOFA_ROOT "${SOFA_MODULE_DIR}/../.." ABSOLUTE)
endif()

if (NOT SOFA_FIND_COMPONENTS)
    set(SOFA_FIND_COMPONENTS
            # SofaCommon
            SofaSimpleFem SofaRigid SofaDeformable SofaObjectInteraction SofaMeshCollision SofaEngine SofaExplicitOdeSolver SofaImplicitOdeSolver SofaLoader

            # SofaBase
            SofaBaseLinearSolver SofaEigen2Solver SofaBaseTopology SofaBaseCollision SofaBaseMechanics SofaBaseVisual SofaBaseUtils

            # SofaMisc
            SofaMiscExtra SofaMiscEngine SofaMiscFem SofaMiscForceField SofaMiscMapping SofaMiscSolver SofaMiscTopology

            # SofaAdvanced
            SofaNonUniformFem

            # SofaGeneral
            SofaBoundaryCondition SofaGeneralMeshCollision SofaGeneralVisual SofaGraphComponent SofaGeneralAnimationLoop
            SofaGeneralDeformable SofaGeneralEngine SofaGeneralExplicitOdeSolver SofaGeneralImplicitOdeSolver
            SofaGeneralLinearSolver SofaGeneralRigid SofaGeneralObjectInteraction SofaGeneralSimpleFem SofaGeneralTopology
            SofaTopologyMapping SofaUserInteraction SofaConstraint SofaGeneralLoader

            # Not yet pluginized
            SofaFramework SofaGui SofaSimulation)
endif()

set(_COMPONENT_FOUND )
foreach(component ${SOFA_FIND_COMPONENTS})
    find_package(${component} CONFIG QUIET HINTS ${SOFA_MODULE_DIR})
    list(APPEND _COMPONENT_FOUND "${component}_FOUND")
endforeach()

# Handle the QUIETLY and REQUIRED arguments and set SOFA_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SOFA DEFAULT_MSG ${_COMPONENT_FOUND})
unset(_COMPONENT_FOUND)
MARK_AS_ADVANCED(SOFA_MODULE_DIR SOFA_ROOT)

