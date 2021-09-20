#!/usr/bin/python3

import Sofa

newton_iterations = 20
cg_iterations = 1000
increments=100

cg_precond = 'Diagonal'

poissonRatio = 0.38
youngModulus = 1000
mu = youngModulus / (2.0 * (1.0 + poissonRatio))
l = youngModulus * poissonRatio / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio))

def createScene(root):

    root.addObject('APIVersion', level='21.06')
    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaCaribou')
    root.gravity= [0, 0, 0 ]
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    root.addObject('EulerImplicitSolver', rayleighStiffness="0",  rayleighMass="0.1", vdamping="3")
    root.addObject('CGLinearSolver',  threshold="0.000000001", tolerance="0.0000000001", iterations="100", printLog="false")

    """ hexa = beam_rigid.addChild("hexa")
    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, -1.5, -1.50 + (10)], max=[20, 1.5, 1.5 + (10)], n=[10, 2, 2])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('UniformMass', totalMass="20", showAxisSizeFactor="1")
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')
    hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices="0 10 20 30")
    hexa.addObject('BarycentricMapping', template="Vec3d,Vec3d",  name="mapping", input="@./hexa_mo", output="@../beam_vec/beam_mo") """

    beam_rigid = root.addChild("beam_rigid")
    beam_rigid.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', showObject=True, showObjectScale=3, position=[0, 0, 10, 0, 0, 0, 1, 5, 0, 10, 0, 0, 0, 1,  10, 0, 10, 0, 0, 0, 1, 15, 0, 10, 0, 0, 0, 1])#  15 0 0 0 0 0 1  -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
    beam_rigid.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    beam_rigid.addObject('MeshTopology', name='beam_topo', lines="0 1  1 2  2 3 ")
    beam_rigid.addObject('FixedConstraint', name="FixedConstraint", indices="0 ")
    beam_rigid.addObject('BeamFEMForceField', name="FEM", radius="0.5", youngModulus="100000", poissonRatio="0.40")

    beam_vec = beam_rigid.addChild('beam_vec')
    beam_vec.addObject('MechanicalObject', name="beam_mo", template='Vec3d', position=[0, 0, 10, 5, 0, 10, 10, 0, 10, 15, 0, 10])
    beam_vec.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    beam_vec.addObject('SphereCollisionModel', radius='2') 
    #beam_vec.addObject('ConstantForceField', indices="3", force="0 10 0")
    beam_vec.addObject('IdentityMapping')

    target = root.addChild("target")
    target.addObject('MechanicalObject', name="target", template='Vec3d', position=[0, 0, 10, 5, 0, 10, 10, 0, 10, 15, 5, 10])
    target.addObject('SphereCollisionModel', color="blue", radius='1')
    target.addObject('FixedConstraint', name="FixedConstraint", indices="0 1 2 3")
    target.addObject('StiffSpringForceField', template="Vec3", spring="3 3 100 5 0", object1="@../beam_rigid/beam_vec", object2="@./")   
