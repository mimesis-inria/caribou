.. _sphere_iso_doc:
.. role:: important
.. role:: warning

<SphereIsoSurface />
====================

Implicit surface of a sphere to be used with level-set compatible components such as the :ref:`fictitious_grid_doc`.


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
      - Radius of the sphere.
    * - center
      - list
      - [0, 0, 0]
      - Coordinates at the center of the sphere.

Quick example
*************
.. content-tabs::

    .. tab-container:: tab1
        :title: XML

        .. code-block:: xml

            <Node>
                <SphereIsoSurface radius="5" center="0 0 0" />
                <FictitiousGrid name="grid" template="Vec3d" n="4 4 4" min="-5 -5 -5" max="5 5 5" draw_boundary_cells="1" printLog="1" />
            </Node>

    .. tab-container:: tab2
        :title: Python

        .. code-block:: python

            node.addObject('SphereIsoSurface', radius=5, center=[0, 0, 0])
            node.addObject('FictitiousGrid',
                           template='Vec3d',
                           name='grid',
                           n=[4, 4, 4],
                           min=[-5, -5, -5],
                           max=[+5, +5, +5],
                           printLog=True,
                           draw_boundary_cells=True,
                           )

Available python bindings
*************************

None at the moment.
