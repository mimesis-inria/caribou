#!/usr/bin/python3

import Sofa
import SofaCaribou

def addBeam(root, use_caribou, c):
    beam = root.addChild('Beam')
    beam.addObject('EulerImplicitSolver', rayleighStiffness=0, rayleighMass=0)
    beam.addObject('CGLinearSolver', iterations=2500, tolerance=1e-9, threshold=1e-9)

    beam.addObject('RegularGridTopology', name='mesh', min=[c[0]-7.5, c[1]-7.5, c[2]], max=[c[0]+7.5, c[1]+7.5, c[2]+80], n=[5, 5, 21])
    beam.addObject('MechanicalObject', src='@mesh')

    beam.addObject('HexahedronSetTopologyContainer', name='hexahedron_container', src='@mesh')
    beam.addObject('HexahedronSetGeometryAlgorithms')

    # beam.addObject('MeshMatrixMass', massDensity=2.3)

    if use_caribou:
        beam.addObject('HexahedronElasticForce', poissonRatio=0, youngModulus=1e5, linearStrain=True, corotated=True, integration_method="Regular", topology_container='@hexahedron_container', printLog=True)
    else:
        beam.addObject('HexahedronFEMForceField', poissonRatio=0, youngModulus=1e5, method="polar", topology='@hexahedron_container', printLog=False)

    beam.addObject('BoxROI', name="fixed_roi", box=[c[0]-7.5, c[1]-7.5, c[2]-0.1, c[0]+7.5, c[1]+7.5, c[2]+0.1], drawBoxes=True)
    beam.addObject('FixedConstraint', indices="@fixed_roi.indices")

    beam.addObject('BoxROI', name='right_roi',
                   strict=True,
                   box=[c[0]-7.5, c[1]-7.5, c[2]+80-0.1, c[0]+7.5, c[1]+7.5, c[2]+80+0.1],
                   drawBoxes=True)
    beam.addObject('QuadSetTopologyContainer', name='quads_container', quads='@right_roi.quadInROI')
    beam.addObject('TriangleSetTopologyContainer', name='triangles_container')
    beam.addObject('TriangleSetTopologyModifier')
    beam.addObject('Quad2TriangleTopologicalMapping', input='@quads_container', output='@triangles_container')
    beam.addObject('TractionForce', traction=[0, 0, 30000], slope=1/50., triangles='@triangles_container.triangles')

def createScene(root):
    root.dt = 1
    root.addObject('APIVersion', level='17.06')
    root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")

    addBeam(root, True, [0, 0, 0])
    addBeam(root, False, [0, 45, 0])


if __name__ == "__main__":
    import SofaRuntime
    SofaRuntime.importPlugin('SofaComponentAll')
    root = Sofa.Core.Node()
    createScene(root)
    Sofa.Simulation.init(root)
    Sofa.Simulation.animate(root, 1)



