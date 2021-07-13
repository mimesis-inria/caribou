################################################################################
#
# \file      cmake/FindMKL.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Math Kernel Library from Intel
# \date      Thu 26 Jan 2017 02:05:50 PM MST
#
################################################################################

# Find the Math Kernel Library from Intel
#
#  INPUT:
#  MKL_THREADING_VENDOR - GNU, INTEL, SEQUENTIAL or AUTO (default)
#  MKL_STATIC - Find the static libraries instead of the shared ones.
#
#  OUTPUT:
#  MKL_FOUND - System has MKL
#  MKL_INCLUDE_DIRS - MKL include files directories
#  MKL_LIBRARIES - The MKL libraries
#  MKL_FLAGS - The MKL ABI flags
#  MKL_INTERFACE_LIBRARY - MKL interface library
#  MKL_THREAD_LAYER_LIBRARY - MKL threading/sequential layer library
#  MKL_CORE_LIBRARY - MKL core library
#
#  The environment variables MKLROOT and INTEL are used to find the library.
#  Everything else is ignored.
#
#  Example usage:
#
#  set(MKL_STATIC ON)
#  set(MKL_THREADING_VENDOR GNU)
#  find_package(MKL)
#  if(MKL_FOUND)
#    target_link_libraries(${PROJECT_NAME} PRIVATE ${MKL_LIBRARIES})
#    target_compile_options(${PROJECT_NAME} PRIVATE ${MKL_FLAGS})
#  endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT APPLE)
    set(INT_LIB "mkl_gf_lp64")
else()
    set(INT_LIB "mkl_intel_lp64")
endif()

set(COR_LIB "mkl_core")
set(THR_LIB "mkl_sequential")

set(MKL_OPENMP_LIBS)
if (NOT MKL_THREADING_VENDOR OR MKL_THREADING_VENDOR STREQUAL "AUTO")
    find_package(OpenMP QUIET)
    if (OPENMP_FOUND)
        foreach(openmp_lib ${OpenMP_CXX_LIB_NAMES})
            if (openmp_lib MATCHES "^gomp")
                set(MKL_THREADING_VENDOR "GNU")
            elseif(openmp_lib MATCHES "^iomp")
                set(MKL_THREADING_VENDOR "INTEL")
            endif()
        endforeach()
        set(MKL_OpenMP_LIBS ${OpenMP_CXX_LIBRARIES})
    endif()
endif()

if (MKL_THREADING_VENDOR)
    if (MKL_THREADING_VENDOR STREQUAL "GNU")
        set(THR_LIB "mkl_gnu_thread")
    elseif (MKL_THREADING_VENDOR STREQUAL "INTEL")
        set(THR_LIB "mkl_intel_thread")
    endif()
endif()

if(MKL_STATIC)
    set(INT_LIB "lib${INT_LIB}.a")
    set(COR_LIB "lib${COR_LIB}.a")
    set(THR_LIB "lib${THR_LIB}.a")
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

find_path(MKL_INCLUDE_DIR
    NAME mkl.h
    HINTS
        $ENV{MKLROOT}/include
        $ENV{HOME}/intel/mkl/include
        /opt/intel/mkl/include
        /opt/intel/oneapi/mkl/latest/include
        /usr/include/mkl)

if (MKL_INCLUDE_DIR)
    set(MKL_ROOT "${MKL_INCLUDE_DIR}/..")
endif()

# This is where the different MKL libraries will be searched for
list( APPEND _MKL_LIBRARY_SEARCH_DIRECTORIES
        ${MKL_ROOT}/lib
        ${MKL_ROOT}/lib/intel64
        ${MKL_ROOT}/lib/intel64_lin
        $ENV{MKLROOT}/lib
        $ENV{MKLROOT}/lib/intel64
        $ENV{MKLROOT}/lib/intel64_lin
        $ENV{HOME}/intel/mkl/lib
        $ENV{HOME}/intel/mkl/lib/intel64
        $ENV{HOME}/intel/mkl/lib/intel64_lin
        /usr/lib/x86_64-linux-gnu
)

find_library(MKL_INTERFACE_LIBRARY
    NAME ${INT_LIB}
    PATHS
        ${_MKL_LIBRARY_SEARCH_DIRECTORIES})

find_library(MKL_THREAD_LAYER_LIBRARY
    NAME ${THR_LIB}
    PATHS
        ${_MKL_LIBRARY_SEARCH_DIRECTORIES})

find_library(MKL_CORE_LIBRARY
    NAME ${COR_LIB}
    PATHS
        ${_MKL_LIBRARY_SEARCH_DIRECTORIES})

set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
if (APPLE)
    set(MKL_LIBRARIES "${MKL_INTERFACE_LIBRARY};${MKL_THREAD_LAYER_LIBRARY};${MKL_CORE_LIBRARY};${MKL_OpenMP_LIBS}")
else()
    set(MKL_LIBRARIES "-Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_THREAD_LAYER_LIBRARY} ${MKL_CORE_LIBRARY} -Wl,--end-group" ${MKL_OpenMP_LIBS})
endif()
set(MKL_LINKER_FLAGS "")
if (MKL_INCLUDE_DIR AND
    MKL_INTERFACE_LIBRARY AND
    MKL_THREAD_LAYER_LIBRARY AND
    MKL_CORE_LIBRARY)

    if (NOT DEFINED ENV{CRAY_PRGENVPGI} AND
        NOT DEFINED ENV{CRAY_PRGENVGNU} AND
        NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
        NOT DEFINED ENV{CRAY_PRGENVINTEL})
        set(MKL_LINKER_FLAGS "-m64")
    endif()
endif()

list(REMOVE_DUPLICATES MKL_LIBRARIES)

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_INTERFACE_LIBRARY MKL_THREAD_LAYER_LIBRARY MKL_CORE_LIBRARY)

MARK_AS_ADVANCED(MKL_INCLUDE_DIRS MKL_LIBRARIES MKL_INTERFACE_LIBRARY MKL_THREAD_LAYER_LIBRARY MKL_CORE_LIBRARY MKL_LINKER_FLAGS MKL_THREADING_VENDOR)
