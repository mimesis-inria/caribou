#!/usr/bin/python3

import Sofa
import SofaCaribou

class SceneController(Sofa.Core.Controller):
    def __init__(self):
        super().__init__()
        self.ff = None

    def onSimulationInitDoneEvent(self, ev):
        print(f"Initial stiffness condition number is '{self.ff.cond()}'")

    def onAnimateEndEvent(self, ev):
        print(f"Current stiffness condition number is '{self.ff.cond()}'")

def createScene(root):
    root.dt = 1
    controller = root.addObject(SceneController())
    root.addObject('APIVersion', level='17.06')
    root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")

    beam = root.addChild('Beam')
    beam.addObject('EulerImplicitSolver', rayleighStiffness=0, rayleighMass=0)
    beam.addObject('CGLinearSolver', iterations=2500, tolerance=1e-9, threshold=1e-9)

    beam.addObject('RegularGridTopology', name='mesh', min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[5, 5, 21])
    beam.addObject('MechanicalObject', src='@mesh')

    beam.addObject('HexahedronSetTopologyContainer', name='hexahedron_container', src='@mesh')
    beam.addObject('HexahedronSetGeometryAlgorithms')

    beam.addObject('MeshMatrixMass', massDensity=2.3)

    controller.ff = beam.addObject('HexahedronElasticForce', name='ff', poissonRatio=0, youngModulus=1e5, linearStrain=True, corotated=True, topology_container='@hexahedron_container', printLog=True)

    beam.addObject('BoxROI', name="fixed_roi", box=[-7.5, -7.5, -0.9, 7.5, 7.5, 0.1], drawBoxes=True)
    beam.addObject('FixedConstraint', indices="@fixed_roi.indices")

if __name__ == "__main__":
    import SofaRuntime
    SofaRuntime.importPlugin('SofaComponentAll')
    root = Sofa.Core.Node()
    createScene(root)
    Sofa.Simulation.init(root)
    Sofa.Simulation.animate(root, 1)



