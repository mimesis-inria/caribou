#!/usr/bin/python3

import Sofa

newton_iterations = 20
cg_iterations = 1000
increments=100

cg_precond = 'Diagonal'


poissonRatio = 0.38
youngModulus = 3500
mu = youngModulus / (2.0 * (1.0 + poissonRatio))
l = youngModulus * poissonRatio / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio))


def createScene(root):
    root.addObject('APIVersion', level='21.06')
    root.gravity = [0, -9.81, 0]
    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaCaribou')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    indices = "1 4 7 10"

    indices_tab = indices.split()
    val_force = "-100"
    forces_ = " "
    for indice in indices_tab:
        forces_ += " 0 " + val_force + " 0 "


    hexa = root.addChild("hexa")
    hexa.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    hexa.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, -1.5, -1.50], max=[20, 1.5, 1.5], n=[10, 2, 2])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices="0 10 20 30 40")
    #hexa.addObject('ConstantForceField', indices=indices_tab, forces=forces_)


    beam = hexa.addChild("beam")
    beam.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', position="  10 0 0 0 0 0 1 ")#  15 0 0 0 0 0 1  -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
    beam.addObject('UniformMass', totalMass="1", showAxisSizeFactor="2")
    beam.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d')
