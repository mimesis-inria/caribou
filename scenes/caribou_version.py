#!/usr/bin/python3

import Sofa

newton_iterations = 20
cg_iterations = 1000
increments=100

cg_precond = 'Diagonal'

poissonRatio = 0.499
youngModulus = 3500
mu = youngModulus / (2.0 * (1.0 + poissonRatio))
l = youngModulus * poissonRatio / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio))

def createScene(root):

    root.addObject('APIVersion', level='21.06')
    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('RequiredPlugin', name='SofaImplicitOdeSolver')
    root.addObject('RequiredPlugin', name='SofaGeneralSimpleFem')
    root.addObject('RequiredPlugin', name='SofaBoundaryCondition')
    root.addObject('RequiredPlugin', name='SofaCaribou')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    root.addObject('EulerImplicitSolver', rayleighStiffness="0",  rayleighMass="0.1")
    root.addObject('CGLinearSolver',  threshold="0.000000001", tolerance="0.0000000001", iterations="10", printLog="false")

    surx = 80
    sury = 15
    surz = 15

    hexa = root.addChild("hexa")
    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[-1.5, -0.5, 0], max=[3, 1, 1], n=[surx, sury, surz])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')
    hexa.addObject('NeoHookeanMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True, file_nodes="/home/sidaty/Desktop/MeshDifferenceSurfaceIntersection/nodes.txt")
    
    fixed = ''
    forces_indices = ''
    partial_fixed = ""
    force_ = ''
    force = 2
    for i in range(sury*surz):
        fixed += str(surx*i) + " "
        partial_fixed += str(surx*i - 1) + " "
    
    for i in range(sury*surz):
        forces_indices += str(surx*i + surx - 1) + " "
        force_ += str(force) + " 0 0 "
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices=fixed)
    hexa.addObject('PartialFixedConstraint', name="partialfix", indices=partial_fixed, fixedDirections="0 1 1" )

    hexa.addObject('ConstantForceField', indices=forces_indices, force=force_)

    """ beam_rigid = root.addChild("beam_rigid")
    beam_rigid.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', position=[-1.5, 0.25, 0.5, 0, 0, 0, 1,  -1, 0.25, 0.5, 0, 0, 0, 1,  -0.5, 0.25, 0.5, 0, 0, 0, 1,  0, 0.25, 0.5, 0, 0, 0, 1,   0.5, 0.25, 0.5, 0, 0, 0, 1,  1, 0.25, 0.5, 0, 0, 0, 1,   1.5, 0.25, 0.5, 0, 0, 0, 1,   2, 0.25, 0.5, 0, 0, 0, 1,   2.5, 0.25, 0.5, 0, 0, 0, 1,   3, 0.25, 0.5, 0, 0, 0, 1])#  15 0 0 0 0 0 1  -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
    beam_rigid.addObject('UniformMass', totalMass="1", showAxisSizeFactor="0.3")
    beam_rigid.addObject('MeshTopology', name='beam_topo', lines="0 1  1 2  2 3 3 4 4 5 5 6")
    beam_rigid.addObject('PartialFixedConstraint', name="FixedConstraint", indices="0 ")
    beam_rigid.addObject('BeamFEMForceField', name="FEM", radius="0.5", youngModulus="62000", poissonRatio="0.40")
    #beam_rigid.addObject('ConstantForceField', indices="0 1 2 3", force="0 0 0 0 100 0")

    beam_vec = beam_rigid.addChild('beam_vec')
    beam_vec.addObject('MechanicalObject', name="beam_mo", template='Vec3d', position=[-1.5, 0.25, 0.5,  -1, 0.25, 0.5,  -0.5, 0.25, 0.5,  0, 0.25, 0.5,   0.5, 0.25, 0.5,  1, 0.25, 0.5,   1.5, 0.25, 0.5,   2, 0.25, 0.5,   2.5, 0.25, 0.5,   3, 0.25, 0.5])
    beam_vec.addObject('UniformMass', totalMass="0", showAxisSizeFactor="1")
    beam_vec.addObject('SphereCollisionModel', radius='0.1') 
    beam_vec.addObject('IdentityMapping')

    target = hexa.addChild("target")
    target.addObject('MechanicalObject', name="target", template='Vec3d', position=[-1.5, 0.25, 0.5,  -1, 0.25, 0.5,  -0.5, 0.25, 0.5,  0, 0.25, 0.5,   0.5, 0.25, 0.5,  1, 0.25, 0.5,   1.5, 0.25, 0.5,   2, 0.25, 0.5,   2.5, 0.25, 0.5,   3, 0.25, 0.5])
    target.addObject('SphereCollisionModel', color="blue", radius='0.1')
    #target.addObject('FixedConstraint', name="FixedConstraint", indices="0 1 2 3")
    for i in range(10):
        target.addObject('StiffSpringForceField', template="Vec3d", spring=str(i) + " " + str(i) + " 10000000 0 0", object1="@../../beam_rigid/beam_vec", object2="@./")
    target.addObject('BarycentricMapping')    """


    """ hexa = root.addChild("hexa1")
    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, 0, 2], max=[3, 0.5, 2.5], n=[50, 8, 8])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('UniformMass', totalMass="20", showAxisSizeFactor="1")
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')
    hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    fixed = ''
    forces_indices = ''
    force_ = ''
    force = 10
    for i in range(64):
        fixed += str(50*i) + " "
    
    for i in range(64):
        forces_indices += str(50*i + 49) + " "
        force_ += str(force) + " 0 0 "
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices=fixed)
    hexa.addObject('ConstantForceField', indices=forces_indices, force=force_) """