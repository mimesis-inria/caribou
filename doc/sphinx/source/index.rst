Overview
========

The caribou project is aimed at multiphysics computation.
It brings a plugin that complements `SOFA multiphysics framework <https://www.sofa-framework.org>`_.
It also provides generic c++ utilities, and SOFA components such as solvers and forcefields.

The project is composed of two modules:

1. The **Caribou library** brings multiple geometric, linear analysis and topological tools that are designed to be as independent as possible from external projects.

2. The **Sofa caribou library** is built on top of the **caribou library**, but brings new components to the SOFA project as a plugin.

.. raw:: html

    <object data="_static/img/caribou.svg" type="image/svg+xml" width="85%"></object>

.. toctree::
    :maxdepth: 2
    :hidden:

    Building <Building.rst>


.. toctree::
    :caption: Force fields
    :maxdepth: 1
    :titlesonly:
    :hidden:

    HexahedronElasticForce  <Forcefield/HexahedronElasticForce.rst>
    TetrahedronElasticForce <Forcefield/TetrahedronElasticForce.rst>
    HyperelasticForcefield  <Forcefield/HyperelasticForcefield.rst>
    TractionForcefield      <Forcefield/TractionForcefield.rst>


.. toctree::
    :caption: Mappings
    :maxdepth: 1
    :titlesonly:
    :hidden:

    CaribouBarycentricMapping <Mapping/CaribouBarycentricMapping.rst>

.. toctree::
    :caption: Mass
    :maxdepth: 1
    :titlesonly:
    :hidden:

    CaribouMass <Mass/CaribouMass.rst>

.. toctree::
    :caption: Materials
    :maxdepth: 1
    :titlesonly:
    :hidden:

    SaintVenantKirchhoffMaterial <Material/SaintVenantKirchhoffMaterial.rst>
    NeoHookeanMaterial <Material/NeoHookeanMaterial.rst>

.. toctree::
    :caption: Time integrators
    :maxdepth: 1
    :titlesonly:
    :hidden:

    BackwardEulerODESolver <Ode/BackwardEulerODESolver.rst>
    StaticODESolver <Ode/StaticODESolver.rst>
    LegacyStaticODESolver <Ode/LegacyStaticODESolver.rst>

.. toctree::
    :caption: Linear solvers
    :maxdepth: 1
    :titlesonly:
    :hidden:

    ConjugateGradientSolver <Solver/ConjugateGradientSolver.rst>
    LLTSolver  <Solver/LLTSolver.rst>
    LDLTSolver <Solver/LDLTSolver.rst>
    LUSolver   <Solver/LUSolver.rst>

.. toctree::
    :caption: Topology
    :maxdepth: 1
    :titlesonly:
    :hidden:

    FictitiousGrid <Topology/FictitiousGrid.rst>
    CircleIsoSurface <Topology/CircleIsoSurface.rst>
    CylinderIsoSurface <Topology/CylinderIsoSurface.rst>
    SphereIsoSurface <Topology/SphereIsoSurface.rst>