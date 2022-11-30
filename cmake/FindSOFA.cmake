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
    list(APPEND CMAKE_PREFIX_PATH ${SOFA_ROOT})
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

            # SofaSimulation
            SofaSimulationCommon SofaSimulationGraph

            # Not yet pluginized
            SofaFramework SofaGui)
endif()

# Compatibility layer
set(SOFA_VERSION ${SofaFramework_VERSION})
foreach(component ${SOFA_FIND_COMPONENTS})
    if (SOFA_VERSION VERSION_LESS "22.06.99")
        string(REGEX REPLACE "Sofa.Simulation.*" "SofaSimulation" component ${component})
    endif()
    if (SOFA_VERSION VERSION_LESS "20.12.99")
        string(REGEX REPLACE "SofaSimpleFem|SofaRigid|SofaDeformable|SofaObjectInteraction|SofaMeshCollision|SofaEngine|SofaExplicitOdeSolver|SofaImplicitOdeSolver|SofaLoader" "SofaCommon" component ${component})

        string(REGEX REPLACE "SofaBase.+" "SofaBase" component ${component})
        string(REGEX REPLACE "SofaBaseLinearSolver|SofaEigen2Solver" "SofaBase" component ${component})

        string(REGEX REPLACE "SofaMisc.+" "SofaMisc" component ${component})

        string(REPLACE "SofaNonUniformFem" "SofaAdvanced" component ${component})

        string(REGEX REPLACE "SofaGeneral.+" "SofaGeneral" component ${component})
        string(REGEX REPLACE "SofaBoundaryCondition|SofaGraphComponent|SofaTopologyMapping|SofaUserInteraction|SofaConstraint" "SofaGeneral" component ${component})
    endif()
    if (SOFA_VERSION VERSION_LESS "20.12.99")
        string(REGEX REPLACE "SofaSimulationCommon|SofaSimulationGraph" "SofaSimulation" component ${component})
    endif()
    list(APPEND SOFA_FIND_COMPATIBLE_COMPONENTS ${component})
endforeach()
list(REMOVE_DUPLICATES SOFA_FIND_COMPATIBLE_COMPONENTS)
set(SOFA_FIND_COMPONENTS "${SOFA_FIND_COMPATIBLE_COMPONENTS}")

set(_COMPONENT_FOUND )
foreach(component ${SOFA_FIND_COMPONENTS})
    find_package(${component} CONFIG QUIET HINTS ${SOFA_MODULE_DIR})
    list(APPEND _COMPONENT_FOUND "${component}_FOUND")
    if(${component}_FOUND)
        set(SOFA_${component}_FOUND ${component}_FOUND)
    endif()
endforeach()

# Handle the QUIETLY and REQUIRED arguments and set SOFA_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SOFA
    REQUIRED_VARS "${_COMPONENT_FOUND}"
    HANDLE_COMPONENTS
)
unset(_COMPONENT_FOUND)
MARK_AS_ADVANCED(SOFA_MODULE_DIR SOFA_ROOT)

