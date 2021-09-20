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
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    root.addObject('EulerImplicitSolver', rayleighStiffness="0",  rayleighMass="0.1")
    root.addObject('CGLinearSolver',  threshold="0.000000001", tolerance="0.0000000001", iterations="100", printLog="false")
    indices = "1 4 7 10"

    indices_tab = indices.split()
    val_force = "-100"
    forces_ = " "
    for indice in indices_tab:
        forces_ += " 0 " + val_force + " 0 "

    """ for i, method in enumerate(['frame']): """
    hexa = root.addChild("hexa")
    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, -1.5, -1.50 + (10)], max=[20, 1.5, 1.5 + (10)], n=[10, 2, 2])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('UniformMass', totalMass="20", showAxisSizeFactor="1")
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')
    hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices="0 10 20 30")
    #hexa.addObject('ConstantForceField', indices=indices_tab, forces=forces_)


    beam = hexa.addChild("beam")
    beam.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', showObject=True, showObjectScale=3, position=[0, 0, 10, 0, 0, 0, 1, 5, 0, 10, 0, 0, 0, 1,  10, 0, 10, 0, 0, 0, 1, 15, 0, 10, 0, 0, 0, 1])#  15 0 0 0 0 0 1  -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
    beam.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    beam.addObject('MeshTopology', name='beam_topo', lines="0 1  1 2 2 3 ")
    beam.addObject('FixedConstraint', name="FixedConstraint", indices="0 ")
    beam.addObject('BeamFEMForceField', name="FEM", radius="1", youngModulus="0.62e4", poissonRatio="0.40" )
    beam.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d', rotation_extraction_method='frame')

    beam = root.addChild("beam2")
    beam.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', showObject=True, showObjectScale=3, position=[0, 20, 10, 0, 0, 0, 1, 5, 20, 10, 0, 0, 0, 1,  10, 20, 10, 0, 0, 0, 1, 15, 20, 10, 0, 0, 0, 1])#  15 0 0 0 0 0 1  -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
    beam.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    beam.addObject('MeshTopology', name='beam_topo', lines="0 1  1 2 2 3 ")
    beam.addObject('FixedConstraint', name="FixedConstraint", indices="0 ")
    beam.addObject('BeamFEMForceField', name="FEM", radius="1", youngModulus="0.62e4", poissonRatio="0.40" )
    #beam.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d', rotation_extraction_method='frame')