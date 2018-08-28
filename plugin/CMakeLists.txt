cmake_minimum_required(VERSION 3.1)
project(pluginCaribou)

set(LIB_NAME caribou)

set(HEADER_FILES
        Engine/CutGridEngine.h
        Event/IterativeSolverEvent.h
        Forcefield/FEMForcefield.h
        Forcefield/IBMForcefield.h
        Forcefield/PressureForcefield.h
        Solver/NewtonRaphsonSolver.h
        Helper/LinearAlgebra.h
        Helper/Hexahedron.h
        Helper/Tetrahedron.h
        Helper/Triangle.h
)

set(SOURCE_FILES
        InitCaribou.cpp
        Engine/CutGridEngine.cpp
        Event/IterativeSolverEvent.cpp
        Forcefield/FEMForcefield.cpp
        Forcefield/IBMForcefield.cpp
        Forcefield/PressureForcefield.cpp
        Solver/NewtonRaphsonSolver.cpp
)

file(GLOB SCENE_FILES
        "scenes/*.py"
        "scenes/*.scn"
        "scenes/*/*.py"
        "scenes/*/*.scn"
)

find_package(SofaGeneral REQUIRED)
find_package(SofaBase REQUIRED)

OPTION(WITH_PYTHON_BINDING "Compile python bindings" ON)
if (WITH_PYTHON_BINDING)
    find_package(SofaPython REQUIRED)

    set(HEADER_FILES ${HEADER_FILES}
            Python/PythonEventBinder.h
            Python/SofaCaribouBindings.h
            Python/Bindings/PythonEventBinder.h
    )

    set(SOURCE_FILES ${SOURCE_FILES}
            Python/PythonEventBinder.cpp
            Python/SofaCaribouBindings.cpp
            Python/Bindings/PythonEventBinder.cpp
    )

    add_definitions(-DWITH_PYTHON_BINDING)
endif ()

OPTION(WITH_PARDISO_SOLVER "Compile the pardiso solver" OFF)
if (WITH_PARDISO_SOLVER)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set ( PARDISO_LIB "" CACHE FILEPATH "set to Pardiso library downloaded from the web" )
    if(PARDISO_LIB STREQUAL "")
        message(FATAL_ERROR "PARDISO_LIB variable is required.")
    endif()

    set(HEADER_FILES ${HEADER_FILES}
            Solver/PardisoSolver.h
            )

    set(SOURCE_FILES ${SOURCE_FILES}
            Solver/PardisoSolver.cpp
            )
endif()

add_library(${LIB_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${SCENE_FILES})
target_link_libraries(${LIB_NAME} SofaComponentGeneral SofaBaseTopology)

if (WITH_PYTHON_BINDING)
    target_link_libraries(${LIB_NAME} SofaPython)
endif ()

if (WITH_PARDISO_SOLVER)
    target_link_libraries(${LIB_NAME} ${PARDISO_LIB} blas lapack gomp)
endif ()

target_include_directories(${LIB_NAME}
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        "$<INSTALL_INTERFACE:include>")

install(TARGETS ${LIB_NAME}
        EXPORT ${LIB_NAME}Config
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

foreach ( file ${HEADER_FILES} )
    get_filename_component( dir ${file} DIRECTORY )
    STRING(REGEX MATCH "^src\\/(.*)$" dira ${dir})
    set(dira ${CMAKE_MATCH_1})
    install( FILES ${file} DESTINATION include/${LIB_NAME}/${dira} )
endforeach()

install(EXPORT ${LIB_NAME}Config
        FILE ${LIB_NAME}Targets.cmake
        DESTINATION lib/cmake
        )
install(EXPORT ${LIB_NAME}Config DESTINATION lib/cmake)
export(TARGETS ${LIB_NAME} FILE ${LIB_NAME}Config.cmake)