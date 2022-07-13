import Sofa
import SofaCaribou
import tempfile
import os
import meshio
import numpy as np


class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, root):
        root.addObject('DefaultVisualManagerLoop')
        root.addObject('DefaultAnimationLoop')
        root.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels showWireframe")
        root.addObject('RequiredPlugin',
                       pluginName="SofaGeneralEngine SofaMiscMapping SofaCaribou CImgPlugin SofaMiscCollision SofaGeneralSimpleFem SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine SofaGeneralDeformable SofaMiscFem SofaMeshCollision SofaGeneralLoader")

        root.addObject("EulerImplicitSolver", name="cg_odesolver", printLog=False, rayleighStiffness=0.1,
                       rayleighMass=0.1)
        root.addObject("CGLinearSolver", iterations=25, name="linear solver", tolerance=1.0e-9, threshold=1.0e-9)

        artery = root.addChild("artery")
        artery.addObject("MeshObjLoader", name="loader", filename="../meshes/cylinder.obj")
        artery.addObject("MechanicalObject", name="dofs", src="@loader")
        artery.addObject("BoxROI", name="left", box=[-0.1, 0, 0, 0.1, 0.75, 1], drawBoxes=True)
        artery.addObject("BoxROI", name="right", box=[2.9, 0, 0, 3, 0.75, 1], drawBoxes=True)
        self.balloon_indices = artery.addObject("BoxROI", name="balloon", box=[1.5, 0, 0, 2.5, 0.75, 1], drawBoxes=True)
        artery.addObject("FixedConstraint", name="right_fixed", indices="@right.indices")
        artery.addObject("FixedConstraint", name="left_fixed", indices="@left.indices")
        artery.addObject("CaribouTopology", template="Triangle", position="@loader.position",
                         indices="@loader.triangles")
        artery.addObject("SaintVenantKirchhoffMaterial", young_modulus="3000", poisson_ratio="0.3")
        artery.addObject('HyperelasticForcefield', printLog=True)
        #
        # artery = root.addChild("artery")
        # artery.addObject("MeshObjLoader", name="loader", filename="../meshes/cylinder.obj")
        # artery.addObject("MechanicalObject", name="dofs", src="@loader")
        # artery.addObject("BoxROI", name="left", box=[-0.1, 0, 0, 0.1, 0.75, 1], drawBoxes=True)
        # artery.addObject("BoxROI", name="right", box=[2.9, 0, 0, 3, 0.75, 1], drawBoxes=True)
        # self.balloon_indices = artery.addObject("BoxROI", name="balloon", box=[1.5, 0, 0, 2.5, 0.75, 1], drawBoxes=True)
        # artery.addObject("FixedConstraint", name="right_fixed", indices="@right.indices")
        # artery.addObject("FixedConstraint", name="left_fixed", indices="@left.indices")
        # artery.addObject("TriangleSetTopologyContainer", src="@loader", name="Container")
        # artery.addObject("TriangleSetTopologyModifier", name="Modifier")
        # artery.addObject("TriangleSetGeometryAlgorithms", name="GeomAlgo")
        # artery.addObject("TriangularFEMForceField", name="fem", youngModulus=3000, poissonRatio=0.3)
        #
        # balloon_node = artery.addChild("balloon_node")
        # balloon_node.addObject("MechanicalObject")
        # balloon_node.addObject("SubsetMapping", indices="@../balloon.indices")
        #
        # self.artery2 = root.addChild("artery2")
        # self.artery2.addObject("MeshObjLoader", name="loader", filename="../meshes/cylinder.obj")
        # self.artery2.addObject("MeshSubsetEngine", name="subset", inputPosition="@loader.position",
        #                        indices="@artery/balloon.indices")
        # self.artery2.addObject("MechanicalObject", name="dofs", position="@subset.position")
        # self.artery2.addObject("FixedConstraint", fixAll=True)
        #
        # artery.addObject("VectorSpringForceField", template="Vec3d",
        #                  object1="@./balloon_node",
        #                  object2="@../artery2", stiffness=1e8,
        #                  viscosity=0)

        # visu = artery.addChild("visu")
        # visu.addObject("OglModel", name="Visual", color="red")
        # visu.addObject("IdentityMapping", input="@..", output="@Visual")
        return root

    def onSimulationInitDoneEvent(self, event):
        pass

    def onAnimateBeginEvent(self, event):
        pass

    def onAnimateEndEvent(self, event):
        normals = self.artery2.loader.normals.value[self.balloon_indices.indices.value, :]
        self.artery2.dofs.position = self.artery2.dofs.position + normals * 0.001
        pass


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
