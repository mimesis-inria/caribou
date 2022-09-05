import Sofa
import SofaCaribou
import numpy as np
import meshio

mesh = meshio.read("../meshes/liver.obj")
indices = mesh.cells_dict["triangle"]


class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, root):
        root.addObject('DefaultVisualManagerLoop')
        root.addObject('DefaultAnimationLoop')
        root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
        root.addObject('RequiredPlugin',
                       pluginName="SofaCaribou CImgPlugin SofaMiscCollision SofaGeneralSimpleFem SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine SofaGeneralDeformable SofaMiscFem SofaMeshCollision SofaGeneralLoader")

        root.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-30",
                       absolute_correction_tolerance_threshold="1e-30",
                       absolute_residual_tolerance_threshold="1e-5",
                       relative_residual_tolerance_threshold="1e-10",
                       printLog="1")
        root.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        root.gravity = [0, 0, 0]
        liver = root.addChild('liver')
        liver.addObject("MeshGmshLoader", name="loader", filename="../meshes/liver-refined-3.msh")
        liver.addObject('MechanicalObject', name="mo", src="@loader")
        liver.addObject('CaribouTopology', name='topology', template="Tetrahedron", position="@loader.position",
                        indices="@loader.tetrahedra")
        liver.addObject('BoxROI', name="BC", box="-5 -1 -2 -3 6 5", drawBoxes=True)
        liver.addObject("ConstantForceField", indices="@BC.indices", totalForce=[0, -1000, 0])
        # liver.addObject("CaribouTopology", name="surface", template="Triangle", indices=indices)
        # liver.addObject('TractionForcefield', traction=[0, -0.1, 0], slope=1, topology='@surface',
        #                           printLog=True)
        # liver.addObject('UniformMass', totalMass="1")

        # liver.addObject('FEniCS_Material', template="Tetrahedron", young_modulus="30",
        #                 poisson_ratio="0.3", material_name="NeoHookean")
        liver.addObject('FEniCS_Material', template="Tetrahedron", bulk_modulus=1000000, a=1180, b=8, a_f=18.5 * 10e4,
                        b_f=16, a_s=2.5 * 10e4, b_s=11.1, a_fs=2160, b_fs=11.4, material_name="Ogden")

        liver.addObject('HyperelasticForcefield_FEniCS', printLog=True)
        liver.addObject('BoxROI', name="FixedIndices", box="-1 -1 -1 5 6 5", drawBoxes=True)
        liver.addObject('FixedConstraint', name="FixedConstraint", indices="@FixedIndices.indices",
                        showObject=True)

        liver_collision = liver.addChild('collision')
        liver_collision.addObject("MeshObjLoader", name="loader", filename="../meshes/liver.obj")
        liver_collision.addObject("MeshTopology", src="@loader")
        self.dofs = liver_collision.addObject("MechanicalObject", name="dofs", template="Vec3", showObject=True,
                                              showObjectScale=10)
        liver_collision.addObject("SubsetMapping")

        # liver_visu = liver.addChild('visu')
        # liver_visu.addObject("OglModel", name="liver_ogl", src="@../loader", color="red")
        # liver_visu.addObject("CaribouBarycentricMapping")

        return root

    def onSimulationInitDoneEvent(self, event):
        pass

    def onAnimateBeginEvent(self, event):
        pass

    def onAnimateEndEvent(self, event):
        current_positions = np.array(self.dofs.position.value.copy().tolist()[62])
        rest_positions = np.array(self.dofs.rest_position.value.copy().tolist()[62])
        displacement = []
        for current_point, initial_point in zip(current_positions, rest_positions):
            if np.linalg.norm(current_point - initial_point) != 0:
                displacement.append(np.linalg.norm(current_point - initial_point))
        print(np.mean(displacement))


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
