# Caribou
The caribou project is a set of tools and a plugin that complement sofa framework. It provides generic c++ utilities, a set of sofa's components such as solvers and forcefields, and a complete python module.

## Build dependencies
To build the project, sofa framework compiled project must be found. If sofa was build from scratch, the build's installation path must be provided. It is generated in the sofa's project with the "make install" command and the generated
output is usually found at the 'sofa_path'/build/install directory.

Since sofa relies on Qt5, the Qt's library must also be found.

## Building
Set the two following environment variable:
SOFA_ROOT : The sofa's installation path (usually 'sofa_path'/build/install). It sould contain the directory lib/cmake.
Qt5_DIR : The Qt's installation path (usually '~/Qt5/version/os/lib/cmake/Qt5')

And than, launch the compilation process with
cmake
make
make install (optional)
 
