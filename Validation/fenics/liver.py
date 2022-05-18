import Sofa
import SofaCaribou

CONTACT = False


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

        liver = root.addChild('liver')
        liver.addObject("EulerImplicitSolver", name="cg_odesolver", printLog=False, rayleighStiffness=0.1,
                        rayleighMass=0.1)
        liver.addObject("CGLinearSolver", iterations=25, name="linear solver", tolerance=1.0e-9, threshold=1.0e-9)
        liver.addObject("MeshGmshLoader", name="loader", filename="mesh/liver.msh")
        liver.addObject('MechanicalObject', name="mo", src="@loader")
        liver.addObject('CaribouTopology', name='topology', template="Tetrahedron", position="@loader.position",
                        indices="@loader.tetrahedra")
        liver.addObject('UniformMass', totalMass="250")

        # liver.addObject('FEniCS_Material', template="Tetrahedron", young_modulus="3000",
        #                 poisson_ratio="0.3", material_name="NeoHookean")
        liver.addObject('FEniCS_Material', template="Tetrahedron", bulk_modulus=1000000, a=1180, b=8, a_f=18.5 * 10e4,
                        b_f=16, a_s=2.5 * 10e4, b_s=11.1, a_fs=2160, b_fs=11.4, material_name="Ogden")

        liver.addObject('HyperelasticForcefield_FEniCS', printLog=True)

        # liver_visu = liver.addChild('visu')
        # liver_visu.addObject("OglModel", name="liver_ogl", src="@../loader", color="red")
        # liver_visu.addObject("CaribouBarycentricMapping")

        if CONTACT:
            root.gravity = [0, 0, -9.81]
            root.dt = 0.1 / 3
            root.addObject("DefaultPipeline", verbose=0, depth=10, draw=0)
            root.addObject("BruteForceBroadPhase")
            root.addObject("BVHNarrowPhase")
            root.addObject("DefaultContactManager", name="Response", response="default")
            root.addObject("MinProximityIntersection", name="Proximity", alarmDistance=0.75, contactDistance=0.5)

            liver_collision = liver.addChild('collision')
            liver_collision.addObject("MeshTopology", src="@../loader")
            liver_collision.addObject("TriangleCollisionModel")
            liver_collision.addObject("LineCollisionModel")
            liver_collision.addObject("PointCollisionModel")

            floor_node = root.addChild("floor")
            floor_node.addObject("MeshSTLLoader", name="loader", filename="../meshes/drap_refined.stl")
            floor_node.addObject("MeshTopology", src="@loader")
            floor_node.addObject("MechanicalObject", src="@loader")
            floor_node.addObject("TriangleCollisionModel", name="FloorTriangleModel", simulated="0", moving="0")
            floor_node.addObject("LineCollisionModel", name="FloorLineModel", simulated="0", moving="0")
            floor_node.addObject("PointCollisionModel", name="FloorPointModel", simulated="0", moving="0")
            floor_node.addObject("OglModel", name="FloorV", src="@loader", texturename="textures/floor.bmp")

        else:
            root.gravity = [0, -900.81, 0]
            liver.addObject('BoxROI', name="FixedIndices", box="-1 2 -2 2 6 3", drawBoxes=False)
            liver.addObject('FixedConstraint', name="FixedConstraint", indices="@FixedIndices.indices",
                            showObject=False)

        return root

    def onSimulationInitDoneEvent(self, event):
        pass

    def onAnimateBeginEvent(self, event):
        pass

    def onAnimateEndEvent(self, event):
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
