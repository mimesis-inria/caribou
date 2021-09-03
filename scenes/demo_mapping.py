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
    root.gravity = [0, 0, 0]
    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaCaribou')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    beam = root.addChild("beam")
    beam.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    beam.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    beam.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, 0, 0], max=[20, 1, 1], n=[4, 2, 2])
    beam.addObject('MechanicalObject', name="beam_mo", template="Vec3d", src='@./mesh')
    beam.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    beam.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    beam.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    beam.addObject('FixedConstraint', name="FixedConstraint", indices="0 4  8 12 ")
    #beam.addObject('ConstantForceField', indices="4 5 6 7 ", forces="-100 0 0  -100 0 0   -100 0 0  -100 0 0 ")

    point = beam.addChild("point")
    point.addObject('MechanicalObject', name="point_mo", template='Rigid3d', position="  4 0.5 0.5 0 0 1 0  8 0.5 0.5 0 0 1 0   11 0.5 0.5 0 0 1 0  ")# 0.5 0.5 12 0 0 1 0   0.5 0.5 16 0 0 1 0 ")
    point.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    point.addObject("SphereCollisionModel", name="sphere", radius="1")
    #point.addObject('MeshTopology', name='beam_topo', lines="0 1  1 2")
    #point.addObject('BeamFEMForceField', name="FEM", radius="0.1", youngModulus="0.62e4", poissonRatio="0.40" )
    point.addObject('ConstantForceField', indices="0 1 2", forces="0 0 0  0 2 0  0 0 0  0 2 0   0 0 0  0 2 0 ")
    point.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d')
    
    """ beam = root.addChild("beam2")
    beam.gravity = [0, 0, 0]
    beam.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    beam.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    beam.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, 20, 0], max=[1, 21, 20], n=[2, 2, 4])
    beam.addObject('MechanicalObject', name="beam_mo", template="Vec3d", src='@./mesh')
    beam.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    beam.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    beam.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    beam.addObject('FixedConstraint', name="FixedConstraint", indices="1 2 3 0 ") """

    #root.addObject('BarycentricMapping', template="Vec3d,Rigid3d", output="@./point/point_mo", input='@./beam/beam_mo')
    #point.addObject('ConstantForceField', force="0 -10 0")
    """ point = root.addChild("point2")
    point.addObject('MechanicalObject', name="point_mo", template='Rigid3', position="   3 3 8 0 0 1 0    ")
    point.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    point.addObject("SphereCollisionModel", name="sphere", radius="1") """
    
    """ point2 = beam.addChild("point2")
    point2.addObject('MechanicalObject', name="point_mo", template='Rigid3', position=" 0.5 0.5 6 0 0 0 1 ")
    point2.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    point2.addObject("SphereCollisionModel", name="sphere", radius="1")
    point2.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d')

    point2 = beam.addChild("point3")
    point2.addObject('MechanicalObject', name="point_mo", template='Rigid3', position=" 0.5 0.5 8 0 0 0 1 ")
    point2.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    point2.addObject("SphereCollisionModel", name="sphere", radius="1")
    point2.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d') """