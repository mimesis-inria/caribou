
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was CaribouConfig.cmake.in                            ########

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

# OPTIONS
set(CARIBOU_USE_FLOAT "OFF")
set(CARIBOU_BUILD_TESTS "OFF")
set(CARIBOU_WITH_SOFA "OFF")
set(CARIBOU_WITH_PYTHON_3 "ON")
set(CARIBOU_WITH_EIGEN_MKL "ON")
set(CARIBOU_WITH_OPENMP "ON")

# COMPONENTS
set(CARIBOU_COMPONENTS "Algebra;Config;Geometry;Mechanics;Topology;Python")

if (NOT Caribou_FIND_COMPONENTS)
    set(Caribou_FIND_COMPONENTS ${CARIBOU_COMPONENTS})
endif()

foreach(component ${Caribou_FIND_COMPONENTS})
    if (NOT ";${CARIBOU_COMPONENTS};" MATCHES ${component})
        set(Caribou_FOUND False)
        set(Caribou_NOT_FOUND_MESSAGE "Unsupported component: ${component}. Available components are ${CARIBOU_COMPONENTS}")
    else()
        # For requested component, execute its "config" script
        set_and_check(config_file ${CMAKE_CURRENT_LIST_DIR}/${component}Config.cmake)
        include(${config_file})
        set(Caribou_${component}_FOUND True)
    endif()
endforeach()

check_required_components(Caribou)
