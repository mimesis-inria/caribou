# Install script for directory: /home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Python/Config

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0"
         RPATH "$ORIGIN/../../../:$loader_path/../../../")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou" TYPE SHARED_LIBRARY FILES "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0"
         OLD_RPATH "::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "$ORIGIN/../../../:$loader_path/../../../")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so.21.6.0")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibrariesx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so"
         RPATH "$ORIGIN/../../../:$loader_path/../../../")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou" TYPE SHARED_LIBRARY FILES "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so"
         OLD_RPATH "::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "$ORIGIN/../../../:$loader_path/../../../")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou/Caribou.cpython-38-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheadersx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Caribou/Python" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/cmake-build-debug/include/Caribou/Python/Caribou.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3/site-packages/Caribou" TYPE FILE FILES "/home/sidaty/Desktop/external_plugins/caribou/src/Caribou/Python/Config//__init__.py")
endif()

