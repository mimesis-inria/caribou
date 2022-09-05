#!/usr/bin/python3

import Sofa


class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, root):
        root.dt = 1

        root.addObject('RequiredPlugin', pluginName=[
            'SofaSparseSolver', 'SofaTopologyMapping', 'SofaBoundaryCondition', 'SofaEngine', 'SofaImplicitOdeSolver'
        ])
        root.addObject('RequiredPlugin',
                       pluginName="SofaCaribou CImgPlugin SofaMiscCollision SofaGeneralSimpleFem SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine SofaGeneralDeformable SofaMiscFem SofaMeshCollision SofaGeneralLoader")


        root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showForceFields')

        root.addObject('RegularGridTopology', name='grid', min=[-7.5, -7.5, 0], max=[7.5, 7.5, 80], n=[3, 3, 9])

        meca = root.addChild("mechanical")
        meca.addObject('BackwardEulerODESolver', newton_iterations=10, rayleigh_stiffness=0, rayleigh_mass=0, residual_tolerance_threshold=1e-5, pattern_analysis_strategy="ALWAYS", printLog=True)
        meca.addObject('LDLTSolver', backend="Eigen")
        meca.addObject('MechanicalObject', name='mo', src='@../grid')

        # Complete hexa container
        meca.addObject('HexahedronSetTopologyContainer', src='@../grid', name='mechanical_topology')

        # - Mechanics
        meca.addObject('SaintVenantKirchhoffMaterial', young_modulus=15000, poisson_ratio=0.3)
        meca.addObject('HyperelasticForcefield')

        # - Mass
        meca.addObject('HexahedronSetGeometryAlgorithms')
        meca.addObject('DiagonalMass', massDensity=0.2)

        # Fix the left side of the beam
        meca.addObject('BoxROI', name='fixed_roi', quad='@surface_topology.quad', box=[-7.5, -7.5, -0.9, 7.5, 7.5, 0.1])
        meca.addObject('FixedConstraint', indices='@fixed_roi.indices')


def createScene(node):
    node.addObject(ControlFrame(node))

# Choose in your script to activate or not the GUI
USE_GUI = True


def main():
    import SofaRuntime
    import Sofa.Gui
    SofaRuntime.importPlugin("SofaOpenglVisual")
    SofaRuntime.importPlugin("SofaImplicitOdeSolver")
    SofaRuntime.importPlugin("SofaLoader")

    root = Sofa.Core.Node("root")
    createScene(root)
    Sofa.Simulation.init(root)

    if not USE_GUI:
        for iteration in range(10):
            Sofa.Simulation.animate(root, root.dt.value)
    else:
        Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
        Sofa.Gui.GUIManager.createGUI(root, __file__)
        Sofa.Gui.GUIManager.SetDimension(1080, 1080)
        Sofa.Gui.GUIManager.MainLoop(root)
        Sofa.Gui.GUIManager.closeGUI()


# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()