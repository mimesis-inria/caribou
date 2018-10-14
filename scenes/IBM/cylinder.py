import os
import Sofa
from math import pi as PI
from Caribou import Mesh
# if 'SOFAPLUGIN_PATH' in os.environ:
#     for pluginPath in os.environ['SOFAPLUGIN_PATH'].split(':'):
#         Sofa.addPluginRepository(pluginPath)

l = 100.
w = l/7.

mesh   = Mesh.cylinder(center1=(0, 0, 0), center2=(l, 0, 0), radius=w/2., size=20, dimension=2, quads=False)

def createScene(root):
    root.createObject('APIVersion', name=17.12)
    root.createObject('VisualStyle', displayFlags="showBehaviorModels showCollisionModels showWireframe")
    root.createObject('RequiredPlugin', name='CGALPlugin')
    root.gravity = [0, -9.8, 0]
    root.dt=0.001
    # root.createObject('FreeMotionAnimationLoop')

    cylnode = root.createChild('cylinder')
    cylnode.createObject('EulerImplicitSolver', rayleighStiffness="0", rayleighMass="0")
    cylnode.createObject('CGLinearSolver', iterations="2500", tolerance="1e-09", threshold="1e-09")
    # cylnode.createObject('CentralDifference')
    # surface = cylnode.createObject('Mesh',
    #                      name='cylinder_mesh',
    #                      edges=mesh.surface.edges.tolist(),
    #                      triangles=mesh.surface.triangles.tolist(),
    #                      quads=mesh.surface.quads.tolist(),
    #                      position=mesh.vertices.tolist(),
    #                      )
    # cylnode.createObject('TriangleSetTopologyContainer', src='@cylinder_mesh')

    size = 1
    grid = cylnode.createObject('MultilevelSparseGridTopology',
                                n=[int(round(l/size))+1, int(round(w/size))+1, int(round(w/size))+1],
                                min=[0, -w/2., -w/2.],
                                max=[l, w/2., w/2.],
                                vertices=mesh.vertices.tolist(),
                                input_triangles=mesh.surface.triangles.tolist())
    # cylnode.createObject('MechanicalObject', showObject=True, src=grid.getLinkPath())
    # container = cylnode.createObject('HexahedronSetTopologyContainer', src=grid.getLinkPath())
    # cylnode.createObject('HexahedronSetGeometryAlgorithms')

    # engine = cylnode.createObject('CutGridEngine',
    #                               grid_topology=grid.getLinkPath(),
    #                               # surface_topology=surface.getLinkPath(),
    #                               showHexahedrons="ALL",
    #                               showBoundaryColor=[0, 0, 1, 0.5],
    #                               showOutsideColor=[0, 1, 0, 0.5],
    #                               showInsideColor=[1, 0, 0, 0.5],
    #                               showTriangles=True,
    #                               showTrianglesColor=[1, 0, 1, 0.5],
    #                               )
    # cylnode.createObject('IBMForcefield',
    #                      youngModulus=500000,
    #                      poissonRatio=0.3,
    #                      printLog=True,
    #                      # hexahedrons_flags=engine.getLinkPath() + '.hexahedrons_flags'
    #                      )

    # cylnode.createObject('HexahedronFEMForceField',
    #                      youngModulus=500000,
    #                      poissonRatio=0.3,
    #                      method="small",
    #                      printLog=False,
    #                      # hexahedrons_flags=engine.getLinkPath() + '.hexahedrons_flags'
    #                      )

    # cylnode.createObject('BoxROI', name="fixed_roi", box=[-0.1, -w/2.-0.1, -w/2-0.1, 0.1, w/2.+0.1, w/2.+0.1])
    # cylnode.createObject('FixedConstraint', indices="@fixed_roi.indices")

    # cylnode.createObject('DiagonalMass', massDensity="2.5")
    # visualnode = cylnode.createChild('visual')
    # visualnode.createObject('VisualModel',
    #                         src=surface.getLinkPath(),
    #                         color='blue')
    #
    visualnode = cylnode.createChild('visual2')
    visualnode.createObject('VisualModel',
                            position=mesh.vertices.tolist(),
                            triangles=mesh.surface.triangles.tolist(),
                            color='red')
