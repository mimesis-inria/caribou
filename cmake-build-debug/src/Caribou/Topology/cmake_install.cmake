# Install script for directory: /home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/include/Caribou/Topology/config.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/BarycentricContainer.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/BaseMesh.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/BaseDomain.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/Domain.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology/Grid" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/Grid/Grid.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology/Grid/Internal" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/Grid/Internal/BaseGrid.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology/Grid/Internal" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/Grid/Internal/BaseMultidimensionalGrid.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology/Grid/Internal" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/Grid/Internal/BaseUnidimensionalGrid.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/HashGrid.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Topology" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Topology/Mesh.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou" TYPE FILE FILES
    "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/cmake/Caribou/TopologyConfig.cmake"
    "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/cmake/Caribou/TopologyConfigVersion.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou/TopologyTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou/TopologyTargets.cmake"
         "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/src/Caribou/Topology/CMakeFiles/Export/lib/cmake/Caribou/TopologyTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou/TopologyTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou/TopologyTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/src/Caribou/Topology/CMakeFiles/Export/lib/cmake/Caribou/TopologyTargets.cmake")
endif()

