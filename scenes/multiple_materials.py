#!/usr/bin/python3

import Sofa


def createScene(root):
    root.addObject('APIVersion', level='17.06')

    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaSparseSolver')
    root.addObject('RequiredPlugin', name='SofaBoundaryCondition')
    root.addObject('RequiredPlugin', name='SofaEngine')

    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showForceFields')

    root.addObject('RegularGridTopology', name='grid', min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])

    meca = root.addChild("mechanical")
    meca.addObject('StaticODESolver', newton_iterations=10, printLog=True)
    meca.addObject('LDLTSolver', backend="Pardiso")
    meca.addObject('MechanicalObject', name='mo', src='@../grid')

    # Complete hexa container (only needed for the BC BoxROI later on, a bug from SOFA's side)
    meca.addObject('HexahedronSetTopologyContainer', src='@../grid', name='mechanical_topology')

    #######################################
    # Material 1 (first half of the beam) #
    #######################################

    # - Mechanics
    meca.addObject('BoxROI', name='material_1_box', hexahedra='@../grid.hexahedra', box=[-7.5, -7.5, -0.9, 7.5, 7.5, 40.1])
    meca.addObject('HexahedronSetTopologyContainer', hexahedra='@material_1_box.hexahedraInROI', name='material_1_hexahedrons')
    meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=300000, poisson_ratio=0.499, name='material_1')
    meca.addObject('HyperelasticForcefield', topology="@material_1_hexahedrons", material='@material_1')
    # - Visual
    meca.addObject('QuadSetTopologyContainer', name='material_1_quads')
    meca.addObject('Hexa2QuadTopologicalMapping', input='@material_1_hexahedrons', output='@material_1_quads')
    meca.addObject('QuadSetTopologyModifier')
    meca.addObject('OglModel', position='@mo.position', quads='@material_1_quads.quads', color='red')

    #######################################
    # Material 2 (second half of the beam)
    #######################################

    # - Mechanics
    meca.addObject('BoxROI', name='material_2_box', hexahedra='@../grid.hexahedra', box=[-7.5, -7.5, 30.9, 7.5, 7.5, 80.1])
    meca.addObject('HexahedronSetTopologyContainer', hexahedra='@material_2_box.hexahedraInROI', name='material_2_hexahedrons')
    meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=3000, poisson_ratio=0.3, name='material_2')
    meca.addObject('HyperelasticForcefield', topology="@material_2_hexahedrons", material='@material_2')
    # - Visual
    meca.addObject('QuadSetTopologyContainer', name='material_2_quads')
    meca.addObject('Hexa2QuadTopologicalMapping', input='@material_2_hexahedrons', output='@material_2_quads')
    meca.addObject('QuadSetTopologyModifier')
    meca.addObject('OglModel', position='@mo.position', quads='@material_2_quads.quads', color='blue')

    # Fix the left side of the beam
    meca.addObject('BoxROI', name='fixed_roi', quad='@../grid.quad', box=[-7.5, -7.5, -0.9, 7.5, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    # Apply traction on the right side of the beam
    meca.addObject('BoxROI', name='top_roi', quad='@../grid.quad', box=[-7.5, -7.5, 79.9, 7.5, 7.5, 80.1])
    meca.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/5, topology='@quad_container')
