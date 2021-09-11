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
    root.gravity = [0, -9.81, 0]
    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaCaribou')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    hexa = root.addChild("hexa")
    hexa.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    hexa.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[-10, -2.5, -2.50], max=[10, 2.5, 2.5], n=[8, 4, 4])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices="0 8 16 24 32 40 48 56 64 72 80 88 96 104 112 120")
    #hexa.addObject('ConstantForceField', indices="4 5 6 7 ", forces="-100 0 0  -100 0 0   -100 0 0  -100 0 0 ")

    beam = hexa.addChild("beam")
    beam.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', position=" -10 1.25 1.25 0 0 1 0  0 1.25 1.25 0 0 1 0  10 1.25 1.25 0 0 1 0  ")# 0.5 0.5 12 0 0 1 0   0.5 0.5 16 0 0 1 0 ")
    beam.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    beam.addObject("SphereCollisionModel", name="sphere", radius="1")
    beam.addObject('FixedConstraint', name="FixedConstraint", indices="0")
    beam.addObject('MeshTopology', name='hexa_topo', lines="0 1 1 2")
    beam.addObject('BeamFEMForceField', name="FEM", radius="0.02", youngModulus="0.62e4", poissonRatio="0.40" )
    beam.addObject('ConstantForceField', indices="1 2", forces="  200 0 0 0 0 0   200 0 0 0 0 0")
    beam.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d')