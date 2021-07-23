 .. _neohookean_material_doc:
 .. role:: important

<NeoHookeanMaterial />
======================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::material::NeoHookeanMaterial`

Implementation of a NeoHookean hyperelastic material.


.. list-table::
    :widths: 1 1 1 100
    :header-rows: 1
    :stub-columns: 0

    * - Attribute
      - Format
      - Default
      - Description
    * - young_modulus
      - float
      - 1000
      - Young's modulus of the material.
    * - poisson_ratio
      - float
      - 0.3
      - Poisson's ratio of the material.

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <RegularGridTopology name="grid" min="-7.5 -7.5 0" max="7.5 7.5 80" n="9 9 21" />
                <MechanicalObject src="@grid" />
                <HexahedronSetTopologyContainer name="topology" src="@grid" />
                <NeoHookeanMaterial young_modulus="3000" poisson_ratio="0.49" />
                <HyperelasticForcefield topology="@topology" template="Hexahedron" printLog="1" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject("RegularGridTopology", name="grid", min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])
            node.addObject("MechanicalObject", src="@grid")
            node.addObject("HexahedronSetTopologyContainer", name="topology", src="@grid")
            node.addObject("NeoHookeanMaterial", young_modulus=3000, poisson_ratio=0.49)
            node.addObject("HyperelasticForcefield", topology="@topology", template="Hexahedron", printLog=True)


Available python bindings
*************************

None at the moment.
