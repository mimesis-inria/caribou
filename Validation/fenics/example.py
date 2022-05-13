from turtle import color
import Sofa
import SofaCaribou
import tempfile
import os
import meshio
import numpy as np

ELEMENT_TYPE = "Tetrahedron"
ELEMENT_APPROXIMATION_DEGREE = 2
MATERIAL_MODEL = "NeoHookean"
# TODO improve the manual permutation for matching the redefinition of the hexahedron

if ELEMENT_TYPE == "Tetrahedron" and ELEMENT_APPROXIMATION_DEGREE == 1:
    element_sofa = "Tetrahedron"
    element_fenics = element_sofa
    mesh = meshio.read("../meshes/beam_p1.vtu")
    indices = np.empty(mesh.cells_dict['tri'].shape)
    indices_sofa = mesh.cells_dict['tri3']
    indices_fenics = indices_sofa
elif ELEMENT_TYPE == "Tetrahedron" and ELEMENT_APPROXIMATION_DEGREE == 2:
    element_sofa = "Tetrahedron10"
    element_fenics = element_sofa
    # mesh = meshio.read("../meshes/p2_from_gmsh.msh")
    mesh = meshio.read("../meshes/beam_p2.vtu")
    indices_sofa = mesh.cells_dict['tetra10']
    indices_fenics = indices_sofa[:, [0, 1, 2, 3, 9, 8, 5, 7, 6, 4]]
elif ELEMENT_TYPE == "Hexahedron" and ELEMENT_APPROXIMATION_DEGREE == 1:
    element_sofa = "Hexahedron"
    element_fenics = element_sofa + "_FEniCS"
    mesh = meshio.read("../meshes/beam_q1.vtu")
    indices_sofa = mesh.cells_dict['hexahedron']
    indices_fenics = indices_sofa[:, [4, 5, 0, 1, 7, 6, 3, 2]]
elif ELEMENT_TYPE == "Hexahedron" and ELEMENT_APPROXIMATION_DEGREE == 2:
    element_sofa = "Hexahedron20"
    element_fenics = "Hexahedron_FEniCS20"
    mesh = meshio.read("../meshes/beam_q2.vtu")
    indices_sofa = mesh.cells_dict['hexahedron20']
    indices_fenics = indices_sofa[:, [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]]


else:
    raise ValueError('The element or the approximation degree is not implemented yet.')

if MATERIAL_MODEL == "SaintVenantKirchhoff" or MATERIAL_MODEL == "NeoHookean":
    material = MATERIAL_MODEL
else:
    raise ValueError('The material model is not implemented yet.')


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
        root.gravity = [0, 0, -9.81]
        root.dt = 0.1 / 3

        root.addObject("DefaultPipeline", verbose=0, depth=10, draw=0)
        root.addObject("BruteForceBroadPhase")
        root.addObject("BVHNarrowPhase")
        root.addObject("DefaultContactManager", name="Response", response="default")
        root.addObject("MinProximityIntersection", name="Proximity", alarmDistance=0.75, contactDistance=0.5)

        liver = root.addChild('liver')
        liver.addObject("EulerImplicitSolver", name="cg_odesolver", printLog=False, rayleighStiffness=0.1,
                             rayleighMass=0.1)
        liver.addObject("CGLinearSolver", iterations=25, name="linear solver", tolerance=1.0e-9, threshold=1.0e-9)
        liver.addObject("MeshGmshLoader", name="loader", filename="mesh/liver.msh")
        liver.addObject('MechanicalObject', name="mo", src="@loader")
        liver.addObject('CaribouTopology', name='topology', template="Tetrahedron", position="@loader.position",
                             indices="@loader.tetrahedra")
        liver.addObject('UniformMass', totalMass="250")

        liver.addObject('FEniCS_Material', template="Tetrahedron", young_modulus="3000",
                              poisson_ratio="0.3", material_name="NeoHookean")

        liver.addObject('HyperelasticForcefield_FEniCS', printLog=True)

        liver_visu = liver.addChild('visu')
        liver_visu.addObject("OglModel", name="liver_ogl", src="@../loader", color="red")
        liver_visu.addObject("CaribouBarycentricMapping")

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