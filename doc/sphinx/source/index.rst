Overview
========

The caribou project is aimed at multiphysics computation.
It brings a plugin that complements `SOFA multiphysics framework <https://www.sofa-framework.org>`_.
It also provides generic c++ utilities, and SOFA components such as solvers and forcefields.

The project is composed of two modules:

1. The **Caribou library** brings multiple geometric, linear analysis and topological tools that are designed to be as independent as possible from external projects.

2. The **Sofa caribou library** is built on top of the **caribou library**, but brings new components to the SOFA project as a plugin.

.. raw:: html

    <object data="_static/img/caribou.svg" type="image/svg+xml"></object>

.. toctree::
    :maxdepth: 2
    :hidden:

    installation
    quickstart


.. toctree::
    :caption: Force fields
    :maxdepth: 1
    :titlesonly:
    :hidden:

    HexahedronElasticForce <hexahedron_elastic_force_doc.rst>
    TetrahedronElasticForce <tetrahedron_elastic_force_doc.rst>
    HyperelasticForceField <hyperelastic_forcefield_doc.rst>
    TractionForce <traction_force_doc.rst>


.. toctree::
    :caption: Materials
    :maxdepth: 1
    :titlesonly:
    :hidden:

    SaintVenantKirchhoff <stvk_material_doc.rst>
    NeoHookean <neohookean_material_doc.rst>

.. toctree::
    :caption: ODEs
    :maxdepth: 1
    :titlesonly:
    :hidden:

    StaticODESolver <static_ode_doc.rst>

.. toctree::
    :caption: Linear solvers
    :maxdepth: 1
    :titlesonly:
    :hidden:

    ConjugateGradient <cg_solver_doc.rst>
    SparseLLT <sparse_llt_doc.rst>
    SparseLDLT <sparse_ldlt_doc.rst>
    SparseLU <sparse_lu_doc.rst>