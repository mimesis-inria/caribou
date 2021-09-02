#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Caribou::Caribou.Python.Config" for configuration "Debug"
set_property(TARGET Caribou::Caribou.Python.Config APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(Caribou::Caribou.Python.Config PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_DEBUG "Python::Python"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0"
  IMPORTED_SONAME_DEBUG "Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS Caribou::Caribou.Python.Config )
list(APPEND _IMPORT_CHECK_FILES_FOR_Caribou::Caribou.Python.Config "${_IMPORT_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0" )

# Import target "Caribou::Caribou.Python.Geometry" for configuration "Debug"
set_property(TARGET Caribou::Caribou.Python.Geometry APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(Caribou::Caribou.Python.Geometry PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_DEBUG "Python::Python"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/python3/site-packages/Caribou/Geometry/Geometry.cpython-38-x86_64-linux-gnu.so.21.6.0"
  IMPORTED_SONAME_DEBUG "Geometry.cpython-38-x86_64-linux-gnu.so.21.6.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS Caribou::Caribou.Python.Geometry )
list(APPEND _IMPORT_CHECK_FILES_FOR_Caribou::Caribou.Python.Geometry "${_IMPORT_PREFIX}/lib/python3/site-packages/Caribou/Geometry/Geometry.cpython-38-x86_64-linux-gnu.so.21.6.0" )

# Import target "Caribou::Caribou.Python.Topology" for configuration "Debug"
set_property(TARGET Caribou::Caribou.Python.Topology APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(Caribou::Caribou.Python.Topology PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_DEBUG "Python::Python"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/python3/site-packages/Caribou/Topology/Topology.cpython-38-x86_64-linux-gnu.so.21.6.0"
  IMPORTED_SONAME_DEBUG "Topology.cpython-38-x86_64-linux-gnu.so.21.6.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS Caribou::Caribou.Python.Topology )
list(APPEND _IMPORT_CHECK_FILES_FOR_Caribou::Caribou.Python.Topology "${_IMPORT_PREFIX}/lib/python3/site-packages/Caribou/Topology/Topology.cpython-38-x86_64-linux-gnu.so.21.6.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
