 .. _hexahedron_elastic_force_doc:
 .. role:: important

<HexahedronElasticForce />
==========================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::forcefield::HexahedronElasticForce`

Implementation of a corotational linear elasticity forcefield for hexahedral topologies.

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
    * - integration_method
      - option
      - Regular
      - Integration method used to integrate the stiffness matrix.

            * **Regular**
                  Regular 8 points gauss integration (default).
            * **OnePointGauss**
                  One gauss point integration at the center of the hexahedron

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
                <HexahedronElasticForce topology_container="@topology" youngModulus="3000" poissonRatio="0.49" corotated="1" printLog="1" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject("RegularGridTopology", name="grid", min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])
            node.addObject("MechanicalObject", src="@grid")
            node.addObject("HexahedronSetTopologyContainer", name="topology", src="@grid")
            node.addObject("HexahedronElasticForce", topology_container="@topology", youngModulus=3000, poissonRatio=0.49, corotated=True, printLog=True)


Available python bindings
*************************

.. py:class:: HexahedronElasticForce

    .. py:class:: GaussNode

        :members:
            - **weight** : Gauss node's weight. :class:`numpy.double`
            - **jacobian_determinant** : Gauss node's Jacobian determinant. :class:`numpy.double`
            - **dN_dx** : Gauss node's shape derivatives w.r.t. the current position vector. :class:`numpy.ndarray`
            - **F** : Gauss node's strain tensor. :class:`numpy.ndarray`

    .. py:function:: gauss_nodes_of(hexahedron_id)

        :param hexahedron_id: Index of the hexahedron in the topology container.
        :type hexahedron_id: int
        :return: Reference to the list of Gauss nodes of the element.
        :rtype: list [:class:`GaussNode`]
        :note: No copy involved.

    .. py:function:: stiffness_matrix_of(hexahedron_id)

        :param hexahedron_id: Index of the hexahedron in the topology container.
        :type hexahedron_id: int
        :return: Reference to the elemental 24x24 tangent stiffness matrix
        :rtype: :class:`numpy.ndarray`
        :note: No copy involved.

        Get the elemental 24x24 tangent stiffness matrix of the given hexahedron.

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
