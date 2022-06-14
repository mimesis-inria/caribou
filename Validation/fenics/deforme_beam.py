import Sofa
import SofaCaribou
import meshio
import numpy as np


class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, node): 
        node.name = "root"

        node.dt = 0.01
        node.gravity = [0, 0, -10]
        node.addObject('RequiredPlugin', name='SofaOpenglVisual')
        node.addObject('RequiredPlugin', name='SofaGeneralSimpleFem')
        node.addObject('RequiredPlugin', name='SofaGeneralLoader')
        node.addObject('RequiredPlugin', name='SofaLoader')
        node.addObject('RequiredPlugin', name='SofaGeneralEngine')
        node.addObject('RequiredPlugin', name='SofaBoundaryCondition')
        node.addObject('RequiredPlugin', name='SofaDeformable')
        node.addObject('RequiredPlugin', name='SofaImplicitOdeSolver')
        node.addObject('RequiredPlugin', name='CGALPlugin')
        node.addObject('RequiredPlugin', name='SofaCaribou')
        node.addObject('RequiredPlugin', name='SofaExporter')

        node.addObject('VisualStyle', displayFlags='showVisualModels showWireframe showBehaviorModels showCollisionModels')
        node.addObject('EulerImplicitSolver', rayleighStiffness="0",  rayleighMass="0.1", vdamping="3")
        node.addObject('CGLinearSolver',  threshold="0.000000001", tolerance="0.0000000001", iterations="100", printLog="false")

        beam = node.addChild('beam')
        beam.addObject('MeshSTLLoader', name='LiverSurface', filename="../meshes/beam_p1.stl")
        
        beam.addObject('FictitiousGrid',   
                        template='Vec3',
                        name='integration_grid',
                        surface_positions="@LiverSurface.position", 
                        surface_triangles="@LiverSurface.triangles",
                        n=self.n,
                        printLog=True
                        )
        beam.addObject('StaticODESolver', newton_iterations=self.newton_iterations, correction_tolerance_threshold=1e-8, residual_tolerance_threshold=1e-8, printLog=True)
        beam.addObject('ConjugateGradientSolver', maximum_number_of_iterations=self.cg_iterations, residual_tolerance_threshold=1e-5, preconditioning_method=self.cg_precond, printLog=False)
        
        
        self.dofs = beam.addObject('MechanicalObject', name="dofs", src="@integration_grid")
        beam.addObject('HexahedronSetTopologyContainer', name='topo', src='@integration_grid')
        beam.addObject('SaintVenantKirchhoffMaterial', young_modulus=self.youngModulus, poisson_ratio=self.poissonRatio)
        beam.addObject('FEniCS_Material', template="Tetrahedron", young_modulus="3000",
                              poisson_ratio="0.3", C01=0.7, C10=-0.55, k=0.001, material_name="NeoHookean", path="/home/..")
        beam.addObject('HyperelasticForcefield_FEniCS', printLog=True)

        beam.addObject('FixedConstraint', name="FixedConstraint", indices=" 0 1 2 3 4 5 6")
        
        
        
        visu = beam.addChild('Visu', tags="Visual")
        visu.addObject('OglModel', name="VisualModel", color="blue", src="@../LiverSurface")
        visu.addObject('BarycentricMapping')


        visu.addObject("STLExporter", filename="output", listening="true", printLog='true', exportEveryNumberOfSteps='3')



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
