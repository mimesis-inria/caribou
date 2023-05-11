#!/usr/bin/python3

import Sofa

cg_iterations = 1000
newton_iterations = 20
increments=100
poissonRatio = 0
youngModulus = 3000


def createScene(root):
    root.addObject('RequiredPlugin', pluginName=[
        'Sofa.Component.SceneUtility', # APIVersion
        'Sofa.Component.Constraint.Projective', # FixedConstraint
        'Sofa.Component.Engine.Select', # BoxROI
        'Sofa.Component.LinearSolver.Iterative', # CGLinearSolver
        'Sofa.Component.ODESolver.Backward', # StaticSolver
        'Sofa.Component.StateContainer', # MechanicalObject
        'Sofa.Component.Topology.Container.Dynamic', # HexahedronSetTopologyContainer, QuadSetTopologyContainer
        'Sofa.Component.Topology.Container.Grid', # RegularGridTopology
        'Sofa.Component.Visual' # VisualStyle
    ])
    root.addObject('APIVersion', level='23.06.99')
    root.addObject('DefaultAnimationLoop')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    root.addObject('RegularGridTopology', name='mesh', min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[3, 3, 7])
    # Caribou ODE solver  -  Caribou CG solver  -  No preconditioner
    i = 0
    tx = 20*i
    meca = root.addChild(f"CG_{i}")
    meca.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', preconditioning_method="None", maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-8, printLog=False, verbose=False)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('HexahedronSetTopologyContainer', name='topo', src='@../mesh')
    meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    # meca.addObject('NeoHookeanMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@quad_container', printLog=True)

    # Caribou ODE solver  -  Caribou CG solver  -  Diagonal preconditioner
    i+=1
    tx = 20*i
    meca = root.addChild(f"CG_{i}")
    meca.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', preconditioning_method="Diagonal", maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-8, printLog=False, verbose=False)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('HexahedronSetTopologyContainer', name='topo', src='@../mesh')
    meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    # meca.addObject('NeoHookeanMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@quad_container', printLog=True)

    # SOFA ODE solver  -  Caribou CG solver  -  No preconditioner
    i+=1
    tx = 20*i
    meca = root.addChild(f"CG_{i}")
    meca.addObject('StaticSolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', preconditioning_method="None", maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-8, printLog=False, verbose=False)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('HexahedronSetTopologyContainer', name='topo', src='@../mesh')
    meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    # meca.addObject('NeoHookeanMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@quad_container', printLog=True)

    # SOFA ODE solver  -  Caribou CG solver  -  Diagonal preconditioner
    i+=1
    tx = 20*i
    meca = root.addChild(f"CG_{i}")
    meca.addObject('StaticSolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', preconditioning_method="Diagonal", maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-8, printLog=False, verbose=False)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('HexahedronSetTopologyContainer', name='topo', src='@../mesh')
    meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    # meca.addObject('NeoHookeanMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@quad_container', printLog=True)

    # SOFA ODE solver  -  SOFA CG solver  -  No preconditioner
    i+=1
    tx = 20*i
    meca = root.addChild(f"CG_{i}")
    meca.addObject('StaticSolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('CGLinearSolver', tolerance=1e-8, threshold=1e-25, iterations=cg_iterations)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('HexahedronSetTopologyContainer', name='topo', src='@../mesh')
    meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    # meca.addObject('NeoHookeanMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    meca.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('QuadSetTopologyContainer', name='quad_container', quads='@top_roi.quadInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@quad_container', printLog=True)
