# Caribou
The caribou project is aimed at multiphysics computation. 
It brings a plugin that complements [SOFA multiphysics framework](https://www.sofa-framework.org/). 
It provides generic c++ utilities, a set of sofa components such as solvers and forcefields.

The project is composed of two modules:
1. **Caribou library** brings multiple geometric, linear analysis and topological tools that are designed to 
be as independent as possible from external projects.
2. **Sofa caribou library** is built on top of the **caribou library**, but brings new components to
the SOFA project as a plugin. 

## Build dependencies
To build the **Caribou library**,  a c++17 compliant compiler is required and the external project [Eigen](http://eigen.tuxfamily.org/) has to be installed.

To build the **Sofa caribou library**, SOFA framework binaries and headers must be found. If sofa was build from scratch, 
the build's installation path can be provided to help cmake find those binaries. It is generated in the sofa's project 
with the "make install" command and the generated
output is usually found at the 'sofa_path'/build/install directory.

If your SOFA binaries support Qt5, the Qt's library must also be found.

### Optional dependencies
- The Caribou library also provides python bindings, generated with pybind11.
pybind11 can either be compiled from sources or installed using pip (`sudo python -m pip install pybind11`)
- The unit tests for Caribou depend on the gtest suite (`apt install libgtest-dev` on debian-based linux distros)

## Building

If you are compiling the **Sofa caribou library**, you can set the following environmnent variable
to help cmake find the required binaries:

**SOFA_ROOT** : The sofa's installation path (usually 'sofa_path'/build/install). It should contain the directory lib/cmake.\
**Qt5_DIR** : The Qt's installation path (usually '~/Qt5/version/os/lib/cmake/Qt5')

Launch the compilation process with
```
cmake
make
make install #(optional)
```
