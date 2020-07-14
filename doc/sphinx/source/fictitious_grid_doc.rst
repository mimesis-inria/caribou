.. _fictitious_grid_doc:
.. role:: important
.. role:: note

<FictitiousGrid />
==================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::topology::FictitiousGrid`

Implementation of an advanced fictitious (sparse) grid.

An fictitious grid is a regular grid of hexahedral elements that embed an implicit (iso) or explicit (mesh) surface.
Elements that lie completely outside the embedded surface are ignored, hence the common name of "sparse" grid. This
component allows to retrieve quickly the type (inside, outside or boundary) of a given point location or element. It
also provides a recursive subdivision algorithm of the intersected cells that allow an accurate integration of the
volume inside the embedded surface.

.. raw:: html

    <img width="400px" src="_static/img/fictitious_circle.png" />

:important:`Requires a topology container or an iso-surface.`


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
    * - template
      - str
      - Vec3D
      - Used to specify the world dimension of the grid. (Vec1D, Vec2D or Vec3D)
    * - n
      - [nx, ny, nz]
      -
      - Grid resolution [nx, ny, nz] where nx, ny or nz are the number of nodes in the x,y and z directions.
    * - min
      - [x, y, z]
      -
      - First corner node position of the grid's bounding box. If it is not specified, and an explicit embedded surface
        is given, this will be automatically computed.
    * - max
      - [x, y, z]
      -
      - Second corner node position of the grid's bounding box. If it is not specified, and an explicit embedded surface
        is given, this will be automatically computed.
    * - maximum_number_of_subdivision_levels
      - int
      - 0
      - Number of subdivision levels of the boundary cells (one level split the cell in 4 subcells in 2D, and 8 subcells in 3D).
    * - volume_threshold
      - float
      - 0.0
      - Ignore (tag as outside) every cells having a volume ratio smaller than this threshold.
    * - iso_surface
      - path
      -
      - | Use an implicit surface instead of a tessellated surface. This will be used as a level-set where an iso-value less than zero means the point is inside the boundaries.
        | :note: Cannot be used simultaneously with an explicit surface representation (segments, triangles and quads tesselation)
    * - surface_positions
      - [[x0, y0, z0], ..., [xn, yn, zn]]
      -
      - | Position vector of nodes contained of the explicit embedded surface.
        | :note: Cannot be used simultaneously with an iso_surface.
    * - surface_edges
      - [[e1p1, e1p2], [e2p1, e2p2], ..., [enp1, enp2]]
      -
      - | List of edge's node indices of the explicit embedded surface.
        | :note: Cannot be used simultaneously with an iso_surface.
        | :note: Can only be used in 2D (template="Vec2D")
    * - surface_triangles
      - [[t1p1, t1p2, t1p3], [t2p1, t2p2, t2p3], ..., [tnp1, tnp2, tnp3]]
      -
      - | List of triangle's node indices of the explicit embedded surface.
        | :note: Cannot be used simultaneously with an iso_surface.
        | :note: Can only be used in 3D (template="Vec3D")
    * - surface_quads
      - [[q1p1, q1p2, q1p3, q1p4], [q2p1, q2p2, q2p3, q2p4], ..., [qnp1, qnp2, qnp3, qnp4]]
      -
      - | List of quad's node indices of the explicit embedded surface.
        | :note: Cannot be used simultaneously with an iso_surface.
        | :note: Can only be used in 3D (template="Vec3D")
    * - draw_boundary_cells
      - bool
      - false
      - Draw the cells intersected by the surface boundary.
    * - draw_outside_cells
      - bool
      - false
      - Draw the cells that are outside of the surface boundary.
    * - draw_inside_cells
      - bool
      - false
      - Draw the cells that are inside of the surface boundary.
    * - position
      - [[x0, y0, z0], ..., [xn, yn, zn]]
      -
      - [**OUTPUT**] Position vector of nodes contained in the sparse grid.
    * - quads
      - [[q1p1, q1p2, q1p3, q1p4], [q2p1, q2p2, q2p3, q2p4], ..., [qnp1, qnp2, qnp3, qnp4]]
      -
      - | [**OUTPUT**] List of quads contained in the sparse grid.
        | :note: Only available in 2D (template="Vec2D")
    * - hexahedrons
      - [[h1p1, h1p2, h1p3, h1p4, ..., h1p8], ..., [hnp1, hnp2, hnp3, hnp4, ..., hnp8]]
      -
      - | [**OUTPUT**] List of hexahedrons contained in the sparse grid.
        | :note: Only available in 3D (template="Vec3D")


Quick examples
**************

Using an implicit surface with level-set
----------------------------------------
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <CircleIsoSurface radius="5" center="0 0" />
                <FictitiousGrid name="grid" template="Vec2d" n="4 4" min="-5 -5" max="5 5" maximum_number_of_subdivision_levels="10" draw_boundary_cells="1" draw_outside_cells="1" draw_inside_cells="1" printLog="1" />

                <MechanicalObject template="Vec2d" position="@grid.position" />
                <QuadSetTopologyContainer quads="@grid.quads" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('CircleIsoSurface', radius=5, center=[0, 0])
            node.addObject('FictitiousGrid',
                           template='Vec2d',
                           name='grid',
                           n=[4, 4],
                           min=[-5, -5],
                           max=[+5, +5],
                           maximum_number_of_subdivision_levels=10,
                           printLog=True,
                           draw_boundary_cells=True,
                           draw_outside_cells=True,
                           draw_inside_cells=True
                           )

            node.addObject('MechanicalObject', template='Vec2d', position='@grid.position')
            node.addObject('QuadSetTopologyContainer', quads='@grid.quads')

Using an explicit surface with mesh intersection
------------------------------------------------

.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <MeshVTKLoader name="loader" filename="liver_surface.vtk" />
                <FictitiousGrid name="grid" template="Vec3d" surface_positions="@loader.position" surface_triangles="@loader.triangles" n="20 20 20" maximum_number_of_subdivision_levels="4" draw_boundary_cells="1" printLog="1" />

                <MechanicalObject template="Vec3d" position="@grid.position" />
                <HexahedronSetTopologyContainer hexahedrons="@grid.hexahedrons" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('MeshVTKLoader', name='loader', filename='liver_surface.vtk')
            node.addObject('FictitiousGrid',
                           template='Vec3d',
                           surface_positions='@loader.position',
                           surface_triangles='@loader.triangles'
                           name='grid',
                           n=[20, 20, 20],
                           maximum_number_of_subdivision_levels=4,
                           printLog=True,
                           draw_boundary_cells=True
                           )

            node.addObject('MechanicalObject', template='Vec3d', position='@grid.position')
            node.addObject('HexahedronSetTopologyContainer', hexahedrons='@grid.hexahedrons')

Available python bindings
*************************

.. py:class:: FictitiousGrid

    .. py:function:: number_of_cells()

        :rtype: int

        Get the number of **sparse** cells (inside or on the boundary) in the grid.

    .. py:function:: number_of_nodes()

        :rtype: int

        Get the number of **sparse** nodes (belonging to a **sparse** cell) in the grid.

    .. py:function:: number_of_subdivisions()

        :rtype: int

        Get the number of subdivisions in the grid.

    .. py:function:: cell_volume_ratio_distribution(number_of_decimals=0)

        :param number_of_decimals: Round the volume ratio to the given number of decimals. For example, setting this
                                   value to 2  will generate a distribution of maximum 100 entries (0.00, 0.01, 0.02, ..., 0.99, 1.00).
                                   Setting a value at zero deactivate the rounding of volume ratio.
                                   Default is 0 which means no rounding.
        :type number_of_decimals: int

        :return: A sorted map where the keys are the percentage of volume inside the cell, and the value is a vector
                containing the ids of all cells having this volume percentage.
        :rtype: {:class:`numpy.float`: [int]}

        Compute the distribution of volume ratios of the top level cells of the grid.

        The volume ratio is the ratio of actual volume of a cell over the total volume of the cell.
        Hence, the ratio of a cell outside the boundaries is 0, the ratio of a cell inside is 1,
        and the ratio of boundary cells are between 0 and 1.