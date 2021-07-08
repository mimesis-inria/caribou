 .. _caribou_barycentric_mapping_doc:
 .. role:: important

<CaribouBarycentricMapping />
=============================
.. rst-class:: doxy-label
.. rubric:: Doxygen:
    :cpp:class:`SofaCaribou::mapping::CaribouBarycentricMapping`

Generic barycentric mapping.

The CaribouBarycentricMapping allows to embed nodes into a containing domain. For each embedded nodes, called
mapped nodes, the index of the element (from the container domain) that contains it is stored, paired to the
barycentric coordinates of the node within this element. When paired with a mechanical object (mo), each of
the mo's node positions will automatically follow the parent element that contains it.

Attributes
**********
.. list-table::
    :widths: 1 1 1 100
    :header-rows: 1
    :stub-columns: 0

    * - Attribute
      - Format
      - Default
      - Description
    * - topology
      - path
      - N/A
      - Topology that contains the embedding (parent) elements.

Quick example
*************
Here's an example of a visual model of a cylinder mapped into a rectangular beam. The cylinder is a triangular
mesh, while the rectangular beam is a complete Finite Element solution modelled using a quadratic hexahedral mesh.

.. code-block:: python

    import Sofa, meshio, numpy as np
    from pathlib import Path

    # FE hexahedral mesh
    current_dir = Path(__file__).parent
    beam_q2 = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'beam_q2.vtu').resolve())

    # Mapped surface mesh
    cylinder = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'cylinder_p1.vtu').resolve())

    # Material
    young_modulus = 10000
    poisson_ratio = 0.49

    # Scene creation
    def createScene(root):
        root.addObject('RequiredPlugin', pluginName='SofaCaribou SofaBoundaryCondition SofaEngine SofaOpenglVisual SofaGeneralVisual')
        root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels')
        root.addObject('StaticODESolver', newton_iterations=10, residual_tolerance_threshold=1e-5, pattern_analysis_strategy="BEGINNING_OF_THE_TIME_STEP")
        root.addObject('LDLTSolver', backend="Pardiso")
        root.addChild('mechanics')

        # Mechanical model of the rectangular beam
        root.mechanics.addObject('MechanicalObject', name='mo', position=(mesh.points + p).tolist(), showObject=True, showObjectScale=5)
        root.mechanics.addObject('CaribouTopology', name='volumetric_topology', template=caribou_type, indices=mesh.cells_dict[meshio_type].tolist())
        root.mechanics.addObject('SaintVenantKirchhoffMaterial', young_modulus=young_modulus, poisson_ratio=poisson_ratio)
        root.mechanics.addObject('HyperelasticForcefield')
        root.mechanics.addObject('BoxROI', name='fixed_roi', box=[p[0]-7.5, p[1]-7.5, p[2]-0.9, p[0]+7.5, p[1]+7.5, p[2]+0.1])
        root.mechanics.addObject('FixedConstraint', indices='@fixed_roi.indices')

        # Visual model of the cylinder mapped inside the parent mechanical beam
        root.mechanics.addChild('visual')
        root.mechanics.visual.addObject('CaribouTopology', name='surface_topology', template='Triangle', indices=cylinder.cell_dict['triangle'].tolist(), position=cylinder.points.tolist())
        root.mechanics.visual.addObject('OglModel', name='mo', position='@surface_topology.position', triangles='@surface_topology.indices', color='green')
        root.mechanics.visual.addObject('CaribouBarycentricMapping', topology='../volumetric_topology')

Available python bindings
*************************

None at the moment.
