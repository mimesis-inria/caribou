 .. _traction_force_doc:
 .. role:: important

<TractionForcefield />
======================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::forcefield::TractionForcefield`

Implementation of a traction forcefield for triangle and quad topologies.

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
    * - topology
      - path
      -
      - Path to a topology container (or path to a mesh) that contains the elements on which the traction will be applied.
    * - traction
      - [tx, ty, tz]
      -
      - Tractive force per unit area (if an incremental load is set by the slope parameter, this is the final load
        reached after all increments).
    * - slope
      - float
      - 0
      - Slope of load increment, the resulting tractive force will be p^t = p^{t-1} + p*slope where p is the traction
        force passed as a data and p^t is the traction force applied at time step t. If slope = 0, the traction will be
        constant.
    * - number_of_steps_before_increment
      - int
      - 1
      - Number of time steps to wait before adding an increment. This can be used to simulate Newton-Raphson solving
        process where the time steps are the Newton iterations.
    * - template
      - option
      -
      - The template argument is used to specified the element type on which to compute the traction force.
        By default, the component will try to deduce its element type from the given topology.

            * **Triangle** - Linear triangular elements
            * **Triangle6** - Quadratic triangular elements
            * **Quad** - Linear quadrangle elements
            * **Quad8** - Quadratic quadrangle elements
    * - nodal_forces
      - [[fx, fy, fz], ...]
      -
      - [OUTPUT] Current nodal forces from the applied traction.
    * - total_load
      - float
      - 0
      - [OUTPUT] Accumulated load applied on all the surface area.


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
                <SaintVenantKirchhoffMaterial young_modulus="3000" poisson_ratio="0.49" />
                <HyperelasticForcefield topology="@topology" template="Hexahedron" printLog="1" />

                <BoxROI name="top_roi" box="-7.5 -7.5 79.9 7.5 7.5 80.1" />
                <QuadSetTopologyContainer name="quad_container" quads="@top_roi.quadInROI" />
                <TractionForcefield traction="0 -30 0" slope="0.25" topology="@quad_container" printLog="1" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject("RegularGridTopology", name="grid", min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])
            node.addObject("MechanicalObject", src="@grid")
            node.addObject("HexahedronSetTopologyContainer", name="topology", src="@grid")
            node.addObject("SaintVenantKirchhoffMaterial", young_modulus=3000, poisson_ratio=0.49)
            node.addObject("HyperelasticForcefield", topology="@topology", template="Hexahedron", printLog=True)

            node.addObject('BoxROI', name='top_roi', box=[-7.5, -7.5, 79.9, 7.5, 7.5, 80.1])
            node.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
            node.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@quad_container', printLog=True)


Available python bindings
*************************

None at the moment.