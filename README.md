# Caribou

[![Documentation Status](https://readthedocs.org/projects/caribou/badge/?version=latest)](https://caribou.readthedocs.io/en/latest/?badge=latest)
[![MacOS](https://github.com/jnbrunet/caribou/actions/workflows/macos.yml/badge.svg)](https://github.com/jnbrunet/caribou/actions/workflows/macos.yml)
[![Linux](https://github.com/jnbrunet/caribou/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/jnbrunet/caribou/actions/workflows/ubuntu.yml)

The caribou project is aimed at multiphysics computation. 
It brings a plugin that complements [SOFA multiphysics framework](https://www.sofa-framework.org/). 
It also provides generic c++ utilities, and SOFA components such as solvers and forcefields.

The project is composed of two modules:
1. **Caribou library** brings multiple geometric, linear analysis and topological tools that are designed to 
be as independent as possible from external projects.
2. **Sofa caribou library** is built on top of the **caribou library**, but brings new components to
the SOFA project as a plugin. 

You can read the documentation for building and using Caribou at the following link:
[https://caribou.readthedocs.io](https://caribou.readthedocs.io)
