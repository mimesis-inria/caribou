# Install script for directory: /home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/BaseHexahedron.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/BaseQuad.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/BaseRectangularHexahedron.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/BaseRectangularQuad.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/BaseSegment.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/BaseTetrahedron.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/BaseTriangle.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/Element.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/Hexahedron.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/Quad.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/RectangularHexahedron.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/RectangularQuad.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/Segment.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/Tetrahedron.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Geometry" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Geometry/Triangle.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou" TYPE FILE FILES
    "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/cmake/Caribou/GeometryConfig.cmake"
    "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/cmake/Caribou/GeometryConfigVersion.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou/GeometryTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou/GeometryTargets.cmake"
         "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/src/Caribou/Geometry/CMakeFiles/Export/lib/cmake/Caribou/GeometryTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou/GeometryTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou/GeometryTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Caribou" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/src/Caribou/Geometry/CMakeFiles/Export/lib/cmake/Caribou/GeometryTargets.cmake")
endif()

