.. _cylinder_iso_doc:
.. role:: important
.. role:: warning

<CylinderIsoSurface />
======================

.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::topology::CylinderIsoSurface`

Implicit surface of a circle to be used with level-set compatible components such as the :ref:`fictitious_grid_doc`.


.. list-table::
    :widths: 1 1 1 100
    :header-rows: 1
    :stub-columns: 0

    * - Attribute
      - Format
      - Default
      - Description
    * - radius
      - float
      - 1.0
      - Radius of the cylinder.
    * - length
      - float
      - 5.0
      - length of the cylinder.
    * - center
      - list
      - [0, 0, 0]
      - Coordinates at the center of cylinder.

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <CylinderIsoSurface radius="10" length="200" center="0 0 0" />
                <FictitiousGrid name="grid" template="Vec3d" n="9 9 19" min="-5 -5 -100" max="5 5 100" draw_boundary_cells="1" printLog="1" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('CylinderIsoSurface', radius=10, length=200, center=[0, 0, 0])
            node.addObject('FictitiousGrid',
                           template='Vec3d',
                           name='grid',
                           n=[9, 9, 19],
                           min=[-5, -5, -100],
                           max=[+5, +5, +100],
                           printLog=True,
                           draw_boundary_cells=True,
                           )

Available python bindings
*************************

None at the moment.
