 .. _caribou_mass_doc:
 .. role:: important

<CaribouMass />
===============

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::mass::CaribouMass`

Implementation of a consistent Mass matrix.
The assembly of this mass matrix takes the form of

.. math::
    \boldsymbol{M}_{IK}  = \int_{\Omega_e} \rho_0 N_I N_K d\Omega \boldsymbol{I}

where :math:`I` and :math:`K` are a pair of indices of the element :math:`K` nodes.
Here, :math:`\rho_0` is the mass density as the mass per volume unit (ie
:math:`\frac{m}{v}`) at the undeformed configuration. Finally, :math:`N_I(\boldsymbol{\Psi})`
is the shape function of the :math:`I`th element's node evaluated at
local coordinates :math:`\boldsymbol{\Psi}` relative to the reference (canonical) element.

A diagonal consistent mass matrix is also constructed by scaling down the diagonal
terms in a way that the mass is constant within the element. The procedure is the following:

.. math::
    \boldsymbol{M}_{II}^{\text{diag}}  = s_e M_{II} \boldsymbol{I}  ~ \text{with} ~ M_{II} = \int_{e} \rho_0 N_I^2 d\Omega

With the scaling factor being

.. math::
    s_e  = \frac{M_e}{\sum_I M_{II}} ~\text{, }~ M_e = \int_{e} \rho_0 d\Omega

See the following book for more information:
  Peter Wriggers, Nonlinear finite element methods (2008),
  DOI: 10.1007/978-3-540-71001-1_2

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
    * - lumped
      - bool
      - false
      - Whether or not the mass matrix should be lumped by scaling the diagonal entries such that the mass
        is constant per element. Note that the lumped matrix is always computed. But this parameter
        will determine if it (the lumped) matrix should be used to solve the acceleration (a = M^(-1).f).
    * - density
      - double
      - 1
      - Mass density of the material at the undeformed state formulated as the mass per volume unit,
        ie :math:`\rho_0 = m / v`.
    * - topology
      - path
      -
      - Path to a either a SOFA mesh topology container (such as an `HexahedronSetTopologyContainer` or
        `TetrahedronSetTopologyContainer`) or a CaribouTopology component that contains the elements.
    * - template
      - option
      -
      - The template argument is used to specified the element type on which to compute the mass.
        By default, the component will try to deduce its element type from the given topology.

            * **Tetrahedron** - 4 nodes tetrahedral elements
            * **Tetrahedron10** - 10 nodes tetrahedral elements
            * **Hexahedron** - 8 nodes hexahedral elements
            * **Hexahedron20** - 20 nodes hexahedral elements

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
                <CaribouMass density="2.5" lumped="true" topology="@topology" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject("RegularGridTopology", name="grid", min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])
            node.addObject("MechanicalObject", src="@grid")
            node.addObject("HexahedronSetTopologyContainer", name="topology", src="@grid")
            node.addObject("CaribouMass", density=2.5, lumped=True, topology="@topology")


Available python bindings
*************************

.. py:class:: CaribouMass

    .. py:function:: M()

        :return: Copy of the consistent mass matrix as a compressed column sparse matrix
        :rtype: :class:`scipy.sparse.csc_matrix`
        :note: The mass matrix must have been assembled beforehand. See the assemble_mass_matrix() methods
               to force an assembly.

        Get the consistent mass matrix of a topology as a compressed sparse column major matrix.

    .. py:function:: M_diag()

        :return: Copy of the lumped mass matrix as a compressed column sparse matrix
        :rtype: :class:`scipy.sparse.csc_matrix`
        :note: The mass matrix must have been assembled beforehand. See the assemble_mass_matrix() methods
               to force an assembly.

        The diagonal lumped mass matrix is constructed by scaling down the diagonal terms in a way that the
        mass is constant within the element.

    .. py:function:: assemble(x)

        Assemble the mass matrix M.

        This will force an assembly of the consistent mass matrix. Since the mass matrix is function of
        the the position vector at rest passed as an nx3 array parameter with n the number of nodes.
        If x is omitted, it will use the mechanical state vector "restPosition".

        A copy of the assembled consistent mass matrix M as a column major sparse matrix can be later
        obtained using the method M().