from Caribou import Mesh
import os
import Sofa

if 'SOFAPLUGIN_PATH' in os.environ:
    for pluginPath in os.environ['SOFAPLUGIN_PATH'].split(':'):
        Sofa.addPluginRepository(pluginPath)

# mesh   = Mesh.cylinder(center1=(0, 0, 0), center2=(10, 0, 0), radius=0.5, size=20, dimension=2, quads=True)

def createScene(root):
    root.createObject('APIVersion', name=17.12)
    root.createObject('VisualStyle', displayFlags="showBehaviorModels showCollisionModels showWireframe")
    root.createObject('RequiredPlugin', name='caribou')

    cylnode = root.createChild('cylinder')
    # surface = cylnode.createObject('Mesh',
    #                      name='cylinder_mesh',
    #                      edges=mesh.surface.edges.tolist(),
    #                      triangles=mesh.surface.triangles.tolist(),
    #                      quads=mesh.surface.quads.tolist(),
    #                      position=mesh.vertices.tolist(),
    #                      )
    # cylnode.createObject('TriangleSetTopologyContainer', src='@cylinder_mesh')


    size = 0.25
    grid = cylnode.createObject('RegularGridTopology', n=[int(round(10/size))+1, int(round(1/size))+1, int(round(1/size))+1], min=[0, -0.5, -0.5], max=[10, 0.5, 0.5])
    cylnode.createObject('MechanicalObject', showObject=False, src=grid.getLinkPath())
    container = cylnode.createObject('HexahedronSetTopologyContainer', src=grid.getLinkPath())

    engine = cylnode.createObject('CutGridEngine',
                                  grid_topology=grid.getLinkPath(),
                                  # surface_topology=surface.getLinkPath(),
                                  showHexahedrons="ALL",
                                  showBoundaryColor=[0, 0, 1, 0.5],
                                  showOutsideColor=[0, 1, 0, 0.5],
                                  showInsideColor=[1, 0, 0, 0.5],
                                  showTriangles=True,
                                  showTrianglesColor=[1, 0, 1, 0.5],
                                  )
    cylnode.createObject('IBMForcefield',
                         youngModulus=5000,
                         poissonRatio=0.3,
                         container=container.getLinkPath(),
                         points_flags=engine.getLinkPath() + '.points_flags',
                         hexahedrons_flags=engine.getLinkPath() + '.hexahedrons_flags'
                         )

    # visualnode = cylnode.createChild('visual')
    # visualnode.createObject('VisualModel',
    #                         src=surface.getLinkPath(),
    #                         color='blue')
    #
    # visualnode = cylnode.createChild('visual2')
    # visualnode.createObject('VisualModel',
    #                         position=engine.getLinkPath() + '.triangle_positions',
    #                         triangles=engine.getLinkPath() + '.triangles',
    #                         color='red')
