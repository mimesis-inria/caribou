#!/usr/bin/python3

import Sofa
import meshio
import numpy as np
from pathlib import Path

current_dir = Path(__file__).parent
# FE meshes
beam_p1 = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'beam_p1.vtu').resolve())
beam_p2 = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'beam_p2.vtu').resolve())
beam_q1 = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'beam_q1.vtu').resolve())
beam_q2 = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'beam_q2.vtu').resolve())

meshes = [
    # label, mesh, position, Caribou's element type, meshio's element type
    ("P1", beam_p1, [-10, +50, -100], 'Tetrahedron', 'tetra'),
    ("P2", beam_p2, [-10, -50, -100], 'Tetrahedron10', 'tetra10'),
    ("Q1", beam_q1, [+10, +50, +20], 'Hexahedron', 'hexahedron'),
    ("Q2", beam_q2, [+10, -50, +20], 'Hexahedron20', 'hexahedron20')
]

# Mapped surface
cylinder = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'cylinder_p1.vtu').resolve())
base_indices = cylinder.cells[0].data[np.ma.masked_equal(cylinder.cell_data['gmsh:physical'][0], 1).mask]
top_indices = cylinder.cells[0].data[np.ma.masked_equal(cylinder.cell_data['gmsh:physical'][0], 2).mask]
length_indices = cylinder.cells[0].data[np.ma.masked_equal(cylinder.cell_data['gmsh:physical'][0], 3).mask]

# Material
young_modulus = 10000
poisson_ratio = 0.49


def createScene(root):
    root.addObject('APIVersion', level='21.06')
    root.addObject('RequiredPlugin', pluginName='SofaBoundaryCondition SofaEngine SofaOpenglVisual SofaGeneralVisual')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels')

    root.addObject('StaticODESolver', newton_iterations=10, residual_tolerance_threshold=1e-5, pattern_analysis_strategy="BEGINNING_OF_THE_TIME_STEP")
    root.addObject('LDLTSolver', backend="Pardiso")

    for name, mesh, p, caribou_type, meshio_type in meshes:
        n = root.addChild(name)
        n.addObject('Visual3DText', text=name, position=[p[0], p[1]+20, p[2]-20], color='white', scale=10, depthTest=True)
        n.addObject('MechanicalObject', name='mo', position=(mesh.points + p).tolist(), showObject=True, showObjectScale=5)
        n.addObject('CaribouTopology', name='volumetric_topology', template=caribou_type, indices=mesh.cells_dict[meshio_type].tolist())
        n.addObject('SaintVenantKirchhoffMaterial', young_modulus=young_modulus, poisson_ratio=poisson_ratio)
        n.addObject('HyperelasticForcefield')
        n.addObject('BoxROI', name='fixed_roi', box=[p[0]-7.5, p[1]-7.5, p[2]-0.9, p[0]+7.5, p[1]+7.5, p[2]+0.1])
        n.addObject('FixedConstraint', indices='@fixed_roi.indices')

        n.addChild("top_traction")
        n.top_traction.addObject('CaribouTopology', name='surface_topology', template='Triangle', indices=top_indices.tolist(), position=(cylinder.points + p).tolist())
        n.top_traction.addObject('MechanicalObject', name='mo', position='@surface_topology.position')
        n.top_traction.addObject('TractionForcefield', traction=[0, -500, 0], slope=1/10, topology='@surface_topology')
        n.top_traction.addObject('CaribouBarycentricMapping')

        n.addChild("base_visual")
        n.base_visual.addObject('CaribouTopology', name='surface_topology', template='Triangle', indices=base_indices.tolist(), position=(cylinder.points + p).tolist())
        n.base_visual.addObject('OglModel', name='mo', position='@surface_topology.position', triangles='@surface_topology.indices', color='red')
        n.base_visual.addObject('CaribouBarycentricMapping')

        n.addChild("top_visual")
        n.top_visual.addObject('CaribouTopology', name='surface_topology', template='Triangle', indices=top_indices.tolist(), position=(cylinder.points + p).tolist())
        n.top_visual.addObject('OglModel', name='mo', position='@surface_topology.position', triangles='@surface_topology.indices', color='blue')
        n.top_visual.addObject('CaribouBarycentricMapping')

        n.addChild("length_visual")
        n.length_visual.addObject('CaribouTopology', name='surface_topology', template='Triangle', indices=length_indices.tolist(), position=(cylinder.points + p).tolist())
        n.length_visual.addObject('OglModel', name='mo', position='@surface_topology.position', triangles='@surface_topology.indices', color='green')
        n.length_visual.addObject('CaribouBarycentricMapping')
