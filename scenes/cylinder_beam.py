#!/usr/bin/python3

import Sofa

newton_iterations = 20
cg_iterations = 1000
increments=100

cg_precond = 'Diagonal'


poissonRatio = 0
youngModulus = 3000
mu = youngModulus / (2.0 * (1.0 + poissonRatio))
l = youngModulus * poissonRatio / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio))


def createScene(root):
    root.addObject('APIVersion', level='21.06')
    
    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaCaribou')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    root.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, 0, 0], max=[1, 1, 20], n=[2, 2, 12])

    beam = root.addChild("beam")
    beam.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    beam.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)
    
    beam.addObject('MechanicalObject', name="beam_mo", template="Vec3d", src='@../mesh')
    beam.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@../mesh.hexahedra', position='@../mesh.position')    
    beam.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    beam.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    beam.addObject('FixedConstraint', name="FixedConstraint", indices="1 2 3 0")

    point = beam.addChild("point")
    point.addObject('MechanicalObject', name="point_mo", template='Rigid3', position="0.5 0.5 10 0 0 0 1")
    point.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    point.addObject("SphereCollisionModel", name="sphere", radius="1")
    point.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron')

