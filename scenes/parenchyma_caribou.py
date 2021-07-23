#!/usr/bin/python3

""" import Sofa

newton_iterations = 20
cg_iterations = 1000
increments=100


poissonRatio = 0
youngModulus = 3000

cg_precond = 'Diagonal'

def createScene(root):
    root.addObject('APIVersion', level='21.06')

    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    
    root.addObject('MeshObjLoader', name='LiverSurface', filename='/media/Sidaty/Data/sidaty/external_plugins/caribou/scenes/liver-smooth.obj')

    liver = root.addChild('Liver', gravity='0 -9.81 0')
    liver.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    liver.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)
    liver.addObject('MeshGmshLoader', name='meshLoader', filename='/media/Sidaty/Data/sidaty/external_plugins/caribou/scenes/liver.msh')
    liver.addObject('MechanicalObject', name='dofs', src="@meshLoader")
    liver.addObject('TetrahedronSetTopologyContainer', template='topo', src="@meshLoader")
    liver.addObject('TetrahedronSetTopologyModifier')
    liver.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    liver.addObject('HyperelasticForcefield', template="Tetrahedron", printLog=True)
    liver.addObject('FixedConstraint', name="FixedConstraint", indices="3 40 39")
    
    visu = liver.addChild('Visu', tags="Visual", gravity="0 -9.81 0")
    visu.addObject('OglModel', name="VisualModel", color="red", src="@../../LiverSurface")
    visu.addObject('BarycentricMapping', name="visual mapping", input="@../dofs", output="@VisualModel")

    #visu.addObject('BoxROI', name='top_roi', box=[-7.5, -7.5, 79.9, 7.5, 7.5, 80.1])
    visu.addObject('TriangleSetTopologyContainer', name='triangle_container', triangles='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20')
    visu.addObject('TractionForcefield', traction=[10, -10, 0], slope=1/increments, topology='@triangle_container', printLog=True)

     """


#!/usr/bin/python3

import Sofa

n = [10, 10, 10]
newton_iterations = 20
cg_iterations = 1000
increments=50


poissonRatio = 0
youngModulus = 3000

cg_precond = 'Diagonal'

def createScene(root):
    root.addObject('APIVersion', level='21.06')

    root.addObject('RequiredPlugin', name='SofaComponentAll')
    root.addObject('RequiredPlugin', name='SofaOpenglVisual')
    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    
    root.addObject('MeshObjLoader', name='LiverSurface', filename='/media/Sidaty/Data/sidaty/external_plugins/nv_caribou/scenes/liver-smooth.obj')


    #root.addObject('SphereIsoSurface', radius=radius, center=[0, 0, 0])
    liver = root.addChild('liver', gravity='0 -9.81 0')

    liver.addObject('FictitiousGrid',
                   template='Vec3',
                   name='integration_grid',
                   surface_positions="@LiverSurface.position", 
                   surface_triangles="@LiverSurface.triangles",
                   n=n,
                   printLog=True
                   )
    
    liver.addObject('StaticODESolver', newton_iterations=newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
    liver.addObject('ConjugateGradientSolver', maximum_number_of_iterations=cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=cg_precond, printLog=False)
    
    liver.addObject('MechanicalObject', name="dofs", src="@integration_grid")
    liver.addObject('HexahedronSetTopologyContainer', name='topo', src='@integration_grid')
    liver.addObject('SaintVenantKirchhoffMaterial', young_modulus=youngModulus, poisson_ratio=poissonRatio)
    liver.addObject('HyperelasticForcefield', template="Hexahedron", printLog=True)

    liver.addObject('FixedConstraint', name="FixedConstraint", indices="3 40 39")


    visu = liver.addChild('Visu', tags="Visual", gravity="0 -9.81 0")
    visu.addObject('OglModel', name="VisualModel", color="red", src="@../../LiverSurface")
    visu.addObject('BarycentricMapping', name="visual mapping", input="@../dofs", output="@VisualModel")


    visu.addObject('TriangleSetTopologyContainer', name='triangle_container', triangles='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20')
    visu.addObject('TractionForcefield', traction=[0, -500, 0], slope=1/increments, topology='@triangle_container', printLog=True)

    point = liver.addChild('point')
    point.addObject('MechanicalObject', name="point_mo", template='Vec3d', position="-2 3.5 0 ")
    point.addObject("SphereCollisionModel", name="sphere", color="blue", radius="1")

    point.addObject('BarycentricMapping', name="point mapping", input="@../dofs", output="@./point_mo")
