#!/usr/bin/python3

import Sofa

newton_iterations = 20
cg_iterations = 1000
increments=100

cg_precond = 'Diagonal'


poissonRatio = 0.38
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
    #root.addObject('EulerImplicitSolver', rayleighStiffness="0",  rayleighMass="0.1")
    #root.addObject('CGLinearSolver',  threshold="0.000000001", tolerance="0.0000000001", iterations="100", printLog="false")
    indices = "1 4 7 10"

    indices_tab = indices.split()
    val_force = "-100"
    forces_ = " "
    for indice in indices_tab:
        forces_ += " 0 " + val_force + " 0 "

<<<<<<< HEAD

    hexa = root.addChild("hexa")
    hexa.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    hexa.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, -1.5, -1.50], max=[20, 1.5, 1.5], n=[10, 4, 4])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices="0 9 10 19 20 29 30 39 40 49 50 59 60 69 70 79 80 89 90 99 100 109 110 119 120 129 130 139 140 149 150 159")
    #hexa.addObject('ConstantForceField', indices=indices_tab, forces=forces_)


    beam = hexa.addChild("beam")
    beam.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', position="2.5 0 0 0 0 0 1  5 0 0 0 0 0 1 7.5 0 0 0 0 0 1  10 0 0 0 0 0 1 ")#  15 0 0 0 0 0 1  -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
    beam.addObject('UniformMass', totalMass="1", showAxisSizeFactor="2")
    #beam.addObject('PartialFixedConstraint', indices='0')
    #beam.addObject("SphereCollisionModel", name="sphere", radius="1")
    #beam.addObject('FixedConstraint', name="FixedConstraint", indices="0")
    beam.addObject('MeshTopology', name='hexa_topo', lines="0 1 1 2 2 3")
    beam.addObject('BeamFEMForceField', name="FEM", radius="0.05", radiusInner="0.049",  youngModulus="0.62e4", poissonRatio="0.40")
    #beam.addObject('ConstantForceField', indices="1 2", forces="  200 0 0 0 0 0   200 0 0 0 0 0")
    beam.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d')

    """ hexa = root.addChild("hexa")
    hexa.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    hexa.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[-10, -1.5, -1.50], max=[10, 1.5, 1.5], n=[9, 4, 4])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices="0 9 18 27 36 45 54 63 72 81 90 99 108 117 126 135")
    #hexa.addObject('ConstantForceField', indices=indices_tab, forces=forces_)

    beam = hexa.addChild("beam")
    beam.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', position=" -10 0 0 0 0 0 1  0 0 0 0 0 0 1   10 0 0 0.027 0 0 0 0 0 0 1")# -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
    beam.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    #beam.addObject("SphereCollisionModel", name="sphere", radius="1")
    #beam.addObject('FixedConstraint', name="FixedConstraint", indices="0")
    beam.addObject('MeshTopology', name='hexa_topo', lines="0 1 1 2")
    beam.addObject('BeamFEMForceField', name="FEM", radius="0.05", radiusInner="0.049",  youngModulus="0.62e6", poissonRatio="0.40", useSymmetricAssembly="1")
    #beam.addObject('ConstantForceField', indices="1 2", forces="  200 0 0 0 0 0   200 0 0 0 0 0")
    beam.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d') """

    """ hexa2 = root.addChild("hexa2")
    hexa2.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    hexa2.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    hexa2.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[-1.5, 3, -1.50], max=[1.5, 6, 1.5], n=[4, 4, 4])
    hexa2.addObject('MechanicalObject', name="hexa2_mo", template="Vec3d", src='@./mesh')
    hexa2.addObject('CaribouTopology', name='topo', template='hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    hexa2.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa2.addObject('HyperelasticForcefield', template="hexahedron", printLog=True)
    hexa2.addObject('FixedConstraint', name="FixedConstraint", indices="0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 ")
    hexa2.addObject('ConstantForceField', indices=indices_tab, forces=forces_) """
    
    """ hexa = root.addChild("hexa")
    hexa.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    hexa.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[-10, -2.5, -2.50], max=[10, 2.5, 2.5], n=[8, 4, 4])
    hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
    hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
    hexa.addObject('FixedConstraint', name="FixedConstraint", indices="0 8 16 24 32 40 48 56 64 72 80 88 96 104 112 120")
    hexa.addObject('ConstantForceField', indices=indices_tab, forces=forces_)

    beam = hexa.addChild("beam")
    beam.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', position=" -10 1.25 1.25 0 0 1 0   10 1.25 1.25 0 0 1 0  ")# -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
    beam.addObject('UniformMass', totalMass="1", showAxisSizeFactor="1")
    #beam.addObject("SphereCollisionModel", name="sphere", radius="1")
    #beam.addObject('FixedConstraint', name="FixedConstraint", indices="0")
    beam.addObject('MeshTopology', name='hexa_topo', lines="0 1")
    beam.addObject('BeamFEMForceField', name="FEM", radius="0.01", youngModulus="0.62e6", poissonRatio="0.40" )
    #beam.addObject('ConstantForceField', indices="1 2", forces="  200 0 0 0 0 0   200 0 0 0 0 0")
    beam.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d')

    hexa2 = root.addChild("hexa2")
    hexa2.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    hexa2.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

    hexa2.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[-10, 5, -2.50], max=[10, 10, 2.5], n=[8, 4, 4])
    hexa2.addObject('MechanicalObject', name="hexa2_mo", template="Vec3d", src='@./mesh')
    hexa2.addObject('CaribouTopology', name='topo', template='hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')    
    hexa2.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    hexa2.addObject('HyperelasticForcefield', template="hexahedron", printLog=True)
    hexa2.addObject('FixedConstraint', name="FixedConstraint", indices="0 8 16 24 32 40 48 56 64 72 80 88 96 104 112 120")
    hexa2.addObject('ConstantForceField', indices=indices_tab, forces=forces_)
     """
=======
    for i, method in enumerate(['Frame', 'SVD', 'APD']):
        hexa = root.addChild("hexa")
        hexa.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=False)
        hexa.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)

        hexa.addObject('RegularGridTopology', template="Vec3d", name='mesh', min=[0, -1.5, -1.50 + (i*10)], max=[20, 1.5, 1.5 + (i*10)], n=[10, 2, 2])
        hexa.addObject('MechanicalObject', name="hexa_mo", template="Vec3d", src='@./mesh')
        hexa.addObject('CaribouTopology', name='topo', template='Hexahedron', indices='@./mesh.hexahedra', position='@./mesh.position')
        hexa.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
        hexa.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)
        hexa.addObject('FixedConstraint', name="FixedConstraint", indices="0 10 20 30 40")
        #hexa.addObject('ConstantForceField', indices=indices_tab, forces=forces_)


        beam = hexa.addChild("beam")
        beam.addObject('MechanicalObject', name="beam_mo", template='Rigid3d', showObject=True, showObjectScale=3, position=[10, 0, (i*10), 0, 0, 0, 1,])#  15 0 0 0 0 0 1  -5 1.25 1.25 0 0 1 0  -2.5 1.25 1.25 0 0 1 0   0 1.25 1.25 0 0 1 0  2.5 1.25 1.25 0 0 1 0   5 1.25 1.25 0 0 1 0 ")
        # beam.addObject('UniformMass', totalMass="1", showAxisSizeFactor="2")
        beam.addObject('CaribouBarycentricMapping', topology='@../topo', template='Hexahedron,Rigid3d', rotation_extraction_method=method)
>>>>>>> e59440d270cb56d891f8e333bebd09b9f5bfc69e
