
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was TopologyConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set(CARIBOU_WITH_VTK "OFF")
set(CARIBOU_VTK_MODULES "CommonCore;CommonDataModel;IOLegacy")

find_package(Eigen3 REQUIRED NO_MODULE)

if(CARIBOU_WITH_VTK)
    find_package(VTK COMPONENTS ${CARIBOU_VTK_MODULES} REQUIRED)
    if (VTK_VERSION VERSION_LESS "8.90.0")
        # old system
        include(${VTK_USE_FILE})
    endif()
endif()

if (NOT TARGET Caribou::Topology)
    include("${CMAKE_CURRENT_LIST_DIR}/TopologyTargets.cmake")
endif()

check_required_components(Topology)
