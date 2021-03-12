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

    installation


.. toctree::
    :caption: Force fields
    :maxdepth: 1
    :titlesonly:
    :hidden:

    HexahedronElasticForce <hexahedron_elastic_force_doc.rst>
    TetrahedronElasticForce <tetrahedron_elastic_force_doc.rst>
    HyperelasticForcefield <hyperelastic_forcefield_doc.rst>
    TractionForce <traction_force_doc.rst>


.. toctree::
    :caption: Materials
    :maxdepth: 1
    :titlesonly:
    :hidden:

    SaintVenantKirchhoff <stvk_material_doc.rst>
    NeoHookean <neohookean_material_doc.rst>

.. toctree::
    :caption: ODE
    :maxdepth: 1
    :titlesonly:
    :hidden:

    BackwardEuler <backward_euler_ode_doc.rst>
    StaticODESolver <static_ode_doc.rst>
    LegacyStaticODESolver <legacy_static_ode_doc.rst>

.. toctree::
    :caption: Linear solvers
    :maxdepth: 1
    :titlesonly:
    :hidden:

    ConjugateGradient <cg_solver_doc.rst>
    LLTSolver <sparse_llt_doc.rst>
    LDLTSolver <sparse_ldlt_doc.rst>
    LUSolver <sparse_lu_doc.rst>

.. toctree::
    :caption: Topology
    :maxdepth: 1
    :titlesonly:
    :hidden:

    FictitiousGrid <fictitious_grid_doc.rst>
    CircleIsoSurface <circle_iso_doc.rst>
    CylinderIsoSurface <cylinder_iso_doc.rst>
    SphereIsoSurface <sphere_iso_doc.rst>