 .. _tetrahedron_elastic_force_doc:
 .. role:: important

<TetrahedronElasticForce />
===========================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::forcefield::TetrahedronElasticForce`

Implementation of a corotational linear elasticity forcefield for tetrahedral topologies.

:important:`Requires a mechanical object.`
:important:`Requires a topology container.`


.. list-table::
    :widths: 1 1 1 100
    :header-rows: 1
    :stub-columns: 0

    * - Attribute
      - Format
      - Default
      - Description
    * - printLog
      - bool
      - false
      - Output informative messages at the initialization and during the simulation.
    * - youngModulus
      - float
      - 1000
      - Young's modulus of the material
    * - poissonRatio
      - float
      - 0.3
      - Poisson's ratio of the material
    * - corotated
      - bool
      - true
      - Whether or not to use corotated elements for the strain computation. The rotation is viewed as constant on
        the element and is extracted at its center point.
    * - topology_container
      - path
      -
      - Path to a topology container (or path to a mesh) that contains the hexahedral elements.

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <RegularGridTopology name="grid" min="-7.5 -7.5 0" max="7.5 7.5 80" n="9 9 21">
                <MechanicalObject src="@grid" />
                <HexahedronSetTopologyContainer name="hexahedral_topology" src="@grid" />
                <TetrahedronSetTopologyContainer name="tetrahedral_topology" />
                <TetrahedronSetTopologyModifier />
                <Hexa2TetraTopologicalMapping input="@hexahedral_topology" output="@tetrahedral_topology" />
                <TetrahedronElasticForce topology_container="@tetrahedral_topology" youngModulus="3000" poissonRatio="0.49" corotated="1" printLog="1" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject("RegularGridTopology", name="grid", min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])
            node.addObject("MechanicalObject", src="@grid")
            node.addObject("HexahedronSetTopologyContainer", name="hexahedral_topology", src="@grid")
            node.addObject('TetrahedronSetTopologyContainer', name='tetrahedral_topology')
            node.addObject('TetrahedronSetTopologyModifier')
            node.addObject('Hexa2TetraTopologicalMapping', input='@hexahedral_topology', output='@tetrahedral_topology')
            node.addObject("TetrahedronElasticForce", topology_container="@tetrahedral_topology", youngModulus=3000, poissonRatio=0.49, corotated=True, printLog=True)


Available python bindings
*************************

None at the moment.