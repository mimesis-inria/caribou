#!/usr/bin/python3

import Sofa

newton_iterations = 2
cg_iterations = 100
increments=100

cg_precond = 'None'
# cg_precond = 'Identity'
# cg_precond = 'Diagonal'
# cg_precond = 'LeastSquareDiagonal'
# cg_precond = 'IncompleteCholesky'
# cg_precond = 'IncompleteLU'


poissonRatio = 0
youngModulus = 3000
mu = youngModulus / (2.0 * (1.0 + poissonRatio))
l = youngModulus * poissonRatio / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio))

def createScene(root):
    root.addObject('APIVersion', level='21.06')

    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    root.addObject('RegularGridTopology', name='mesh', min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[9, 9, 21])

    i = 0
    tx = 20*i
    meca = root.addChild("caribou_corotated_tetra")
    meca.addObject('LegacyStaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('TetrahedronSetTopologyContainer', name='topo')
    meca.addObject('TetrahedronSetTopologyModifier')
    meca.addObject('Hexa2TetraTopologicalMapping', input='@../mesh', output='@topo', swapping=True)
    meca.addObject('TetrahedronElasticForce', topology_container='@topo', youngModulus=youngModulus, poissonRatio=poissonRatio, corotated=True, printLog=True)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('TriangleSetTopologyContainer', name='triangle_container', triangles='@top_roi.trianglesInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@triangle_container', printLog=True)


    i += 1
    tx = 20*i
    meca = root.addChild("sofa_corotated_tetra")
    meca.addObject('LegacyStaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    meca.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)
    meca.addObject('MechanicalObject', src='@../mesh', translation=[tx, 0, 0])
    meca.addObject('TetrahedronSetTopologyContainer', name='topo')
    meca.addObject('TetrahedronSetTopologyModifier')
    meca.addObject('Hexa2TetraTopologicalMapping', input='@../mesh', output='@topo', swapping=True)
    meca.addObject('TetrahedronFEMForceField', topology='@topo', youngModulus=youngModulus, poissonRatio=poissonRatio, method='large', printLog=False)

    meca.addObject('BoxROI', name='fixed_roi', box=[-7.5+tx, -7.5, -0.9, 7.5+tx, 7.5, 0.1])
    meca.addObject('FixedConstraint', indices='@fixed_roi.indices')

    meca.addObject('BoxROI', name='top_roi', box=[-7.5+tx, -7.5, 79.9, 7.5+tx, 7.5, 80.1])
    meca.addObject('TriangleSetTopologyContainer', name='triangle_container', triangles='@top_roi.trianglesInROI')
    meca.addObject('TractionForcefield', traction=[0, -30, 0], slope=1/increments, topology='@triangle_container', printLog=True)
