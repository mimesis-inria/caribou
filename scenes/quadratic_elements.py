#!/usr/bin/python3

import Sofa
import meshio
import numpy as np
from pathlib import Path

# Mesh files
current_dir = Path(__file__).parent
beam_p2 = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'beam_p2.vtu').resolve())
beam_q2 = meshio.read((current_dir / '..' / 'Validation' / 'meshes' / 'beam_q2.vtu').resolve())

# Material
young_modulus = 100000
poisson_ratio = 0.49
mass_density  = 2.5


def createScene(root):
    root.addObject('RequiredPlugin', pluginName=[
        'Sofa.Component.SceneUtility', # APIVersion
        'Sofa.Component.Constraint.Projective', # FixedConstraint
        'Sofa.Component.Engine.Select', # BoxROI
        'Sofa.Component.Mass', # UniformMass
        'Sofa.Component.StateContainer', # MechanicalObject
        'Sofa.Component.Visual', # VisualStyle
    ])
    root.addObject('APIVersion', level='23.06.99')
    root.addObject('DefaultAnimationLoop')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    root.addObject('BackwardEulerODESolver', newton_iterations=10, rayleigh_stiffness=0, rayleigh_mass=0, residual_tolerance_threshold=1e-5, pattern_analysis_strategy="ALWAYS", printLog=True)
    root.addObject('LDLTSolver', backend="Eigen")

    root.addChild("tet10")
    root.tet10.addObject('MechanicalObject', name='mo', position=beam_p2.points.tolist(), showObject=True, showObjectScale=5)
    root.tet10.addObject('CaribouTopology', name='volumetric_topology', template='Tetrahedron10', indices=beam_p2.cells_dict['tetra10'].tolist())
    root.tet10.addObject('CaribouMass', density=mass_density, lumped=True, topology='@volumetric_topology')
    root.tet10.addObject('SaintVenantKirchhoffMaterial', young_modulus=young_modulus, poisson_ratio=poisson_ratio)
    root.tet10.addObject('HyperelasticForcefield', printLog=True)
    root.tet10.addObject('BoxROI', name='fixed_roi', box=[-7.5, -7.5, -0.9, 7.5, 7.5, 0.1])
    root.tet10.addObject('FixedConstraint', indices='@fixed_roi.indices')
    root.tet10.addObject('CaribouTopology', name='surface_topology', template='Triangle6', indices=beam_p2.cells_dict['triangle6'][np.array(np.ma.masked_equal(beam_p2.cell_data['gmsh:physical'][0], 2).mask)].tolist())
    root.tet10.addObject('TractionForcefield', traction=[0, 600, 0], slope=1, topology='@surface_topology', printLog=True)

    root.addChild("hex20")
    root.hex20.addObject('MechanicalObject', name='mo', position=beam_q2.points.tolist(), showObject=True, showObjectScale=5, translation=[20, 0, 0])
    root.hex20.addObject('CaribouTopology', name='volumetric_topology', template='Hexahedron20', indices=beam_q2.cells_dict['hexahedron20'].tolist())
    root.hex20.addObject('UniformMass', totalMass=18000*mass_density)  # Volume is 15x15x80 = 18000
    root.hex20.addObject('SaintVenantKirchhoffMaterial', young_modulus=young_modulus, poisson_ratio=poisson_ratio)
    root.hex20.addObject('HyperelasticForcefield', printLog=True)
    root.hex20.addObject('BoxROI', name='fixed_roi', box=[-7.5+20, -7.5, -0.9, 7.5+20, 7.5, 0.1])
    root.hex20.addObject('FixedConstraint', indices='@fixed_roi.indices')
    root.hex20.addObject('CaribouTopology', name='surface_topology', template='Quad8', indices=beam_q2.cells_dict['quad8'][np.array(np.ma.masked_equal(beam_q2.cell_data['gmsh:physical'][0], 2).mask)].tolist())
    root.hex20.addObject('TractionForcefield', traction=[0, 600, 0], slope=1, topology='@surface_topology', printLog=True)
