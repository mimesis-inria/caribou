 .. _hyperelastic_forcefield_doc:
 .. role:: important

<HyperelasticForceField />
==========================

Implementation of an hyperelasticity forcefield for any element type topologies.

:important:`Requires a mechanical object.`
:important:`Requires a topology container.`
:important:`Requires a material.`


.. list-table::
    :widths: 10 10 10 70
    :header-rows: 1
    :stub-columns: 0

    * - Attribute
      - Format
      - Default
      - Description
    * - printLog
      - bool
      - false
      - Whether or not to output informative messages at the initialization and during the simulation.
    * - enable_multithreading
      - bool
      - false
      - Enable the multithreading computation of the stiffness matrix. Only use this if you have a very large number of
        elements, otherwise performance might be worse than single threading. When enabled, use the environment variable
        OMP_NUM_THREADS=N to use N threads.
    * - material
      - path
      -
      - Path to a material component.
    * - topology
      - path
      -
      - Path to a topology container (or path to a mesh) that contains the elements.
    * - draw_scale
      - float
      - 0.85
      - Scaling factor for the drawing of elements (between 0 and 1). The factor allows to shrink the element relative
        to its center point when drawing it.
    * - template
      - option
      -
      - The template argument is used to specified the element type on which to compute the hyperelasticity force.
        By default, the component will try to deduce its element type from the given topology.

            * **Triangle**
            * **Quad**
            * **Tetrahedron**
            * **Hexahedron**

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <RegularGridTopology name="grid" min="-7.5 -7.5 0" max="7.5 7.5 80" n="9 9 21">
                <MechanicalObject src="@grid" />
                <HexahedronSetTopologyContainer name="topology" src="@grid" />
                <SaintVenantKirchhoffMaterial young_modulus="3000" poisson_ratio="0.49" />
                <HyperelasticForcefield topology="@topology" template="Hexahedron" printLog="1" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject("RegularGridTopology", name="grid", min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])
            node.addObject("MechanicalObject", src="@grid")
            node.addObject("HexahedronSetTopologyContainer", name="topology", src="@grid")
            node.addObject("SaintVenantKirchhoffMaterial", young_modulus=3000, poisson_ratio=0.49)
            node.addObject("HyperelasticForcefield", topology="@topology", template="Hexahedron", printLog=True)


Available python bindings
*************************

.. py:class:: HyperelasticForceField

    .. py:function:: K()

        :return: Reference to the forcefield tangent stiffness matrix
        :rtype: :class:`scipy.sparse.csc_matrix`
        :note: No copy involved.

        Get the tangent stiffness matrix of the force field as a compressed sparse column major matrix.

    .. py:function:: cond()

        :return: Condition number of the forcefield's tangent stiffness matrix
        :rtype: :class:`numpy.double`


    .. py:function:: eigenvalues()

        :return: Reference to the eigen values of the forcefield's tangent stiffness matrix.
        :rtype: `list[numpy.double]`
