#!/usr/bin/python3

import Sofa

newton_iterations = 10
cg_iterations = 1000
increments=100
# cg_precond = 'None'
cg_precond = 'Identity'
# cg_precond = 'Diagonal'
# cg_precond = 'LeastSquareDiagonal'
# cg_precond = 'IncompleteCholesky'
# cg_precond = 'IncompleteLU'

poissonRatio = 0
youngModulus = 3000
mu = youngModulus / (2.0 * (1.0 + poissonRatio))
l = youngModulus * poissonRatio / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio))


def createScene(root):
    root.addObject('RequiredPlugin', pluginName=[
        'Sofa.Component.SceneUtility', # APIVersion
        'Sofa.Component.Constraint.Projective', # FixedConstraint
        'Sofa.Component.Engine.Select', # BoxROI
        'Sofa.Component.StateContainer', # MechanicalObject
        'Sofa.Component.Topology.Container.Dynamic', # HexahedronSetTopologyContainer, QuadSetTopologyContainer, TetrahedronSetTopologyContainer, TetrahedronSetTopologyModifier, TriangleSetTopologyContainer
        'Sofa.Component.Topology.Container.Grid', # RegularGridTopology
        'Sofa.Component.Topology.Mapping', # Hexa2TetraTopologicalMapping
        'Sofa.Component.Visual', # VisualStyle
    ])
    root.addObject('APIVersion', level='23.06.99')
    root.addObject('DefaultAnimationLoop')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    root.addObject('RegularGridTopology', name='mesh', min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])

    i = 0
    tx = 20*i
    meca = root.addChild("caribou_StVk_hexa")
    meca.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('HexahedronSetTopologyContainer', name='topo', src='@../mesh')
    # meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('NeoHookeanMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@quad_container', printLog=True)

    i += 1
    tx = 20*i
    meca = root.addChild("caribou_hyperelastic_tetra")
    meca.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('TetrahedronSetTopologyContainer', name='topo')
    meca.addObject('TetrahedronSetTopologyModifier')
    meca.addObject('Hexa2TetraTopologicalMapping', input='@../mesh', output='@topo', swapping=True)
    # meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('NeoHookeanMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('HyperelasticForcefield', template="Tetrahedron", printLog=True)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('TriangleSetTopologyContainer', name='triangle_container', triangles='@top_roi.trianglesInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@triangle_container', printLog=True)
