import Sofa
import SofaCaribou
import meshio
import numpy as np

# ELEMENT_TYPE = "Tetrahedron"
ELEMENT_TYPE = "Hexahedron"
ELEMENT_APPROXIMATION_DEGREE = 2
MATERIAL_MODEL = "NeoHookean"
# MATERIAL_MODEL = "SaintVenantKirchhoff"
TRACTION = [0, -10, 0]
# TODO improve the manual permutation for matching the redefinition of the hexahedron

if ELEMENT_TYPE == "Tetrahedron" and ELEMENT_APPROXIMATION_DEGREE == 1:
    element_sofa = "Tetrahedron"
    element_sonics = element_sofa
    element_surface = "Triangle"
    mesh = meshio.read("../meshes/fenics/cube_tetra_1.msh")
    # mesh = meshio.read("../meshes/fenics/cylinder_tetra_1.msh")
    indices_sofa = mesh.cells_dict['tetra']
    indices_sonics = indices_sofa
    a = [i.size for i in mesh.cell_data["gmsh:physical"]]
    indices = np.cumsum(a)
    left = mesh.cells_dict["triangle"][0:indices[0]]
    right = mesh.cells_dict["triangle"][indices[0]:indices[1]]

    with meshio.xdmf.TimeSeriesReader("../meshes/fenics/cube_tetra_1_solution.xdmf") as reader:
        points, cells = reader.read_points_cells()
        t, point_data, cell_data = reader.read_data(0)

    fenics_rest_position = points
    fenics_rest_position_rearranged = np.empty(fenics_rest_position.shape)
    fenics_position = np.empty(fenics_rest_position.shape)
    # TODO FEniCS and SOFA do not have the same dofs numbering
    from scipy import spatial

    tree = spatial.KDTree(mesh.points)
    for positions, displacement in zip(fenics_rest_position, point_data["f"]):
        fenics_rest_position_rearranged[tree.query(positions)[1]] = positions
        fenics_position[tree.query(positions)[1]] = positions + displacement
    print("initial error:", np.mean(fenics_rest_position_rearranged - mesh.points))

elif ELEMENT_TYPE == "Tetrahedron" and ELEMENT_APPROXIMATION_DEGREE == 2:
    element_sofa = "Tetrahedron10"
    element_sonics = element_sofa
    element_surface = "Triangle6"
    mesh = meshio.read("../meshes/fenics/cube_tetra_2.msh")
    indices_sofa = mesh.cells_dict['tetra10']
    indices_sonics = indices_sofa[:, [0, 1, 2, 3, 9, 8, 5, 7, 6, 4]]
    a = [i.size for i in mesh.cell_data["gmsh:physical"]]
    indices = np.cumsum(a)
    left = mesh.cells_dict["triangle6"][0:indices[0]]
    right = mesh.cells_dict["triangle6"][indices[0]:indices[1]]

    with meshio.xdmf.TimeSeriesReader("../meshes/fenics/cube_tetra_2_solution.xdmf") as reader:
        points, cells = reader.read_points_cells()
        t, point_data, cell_data = reader.read_data(0)

    fenics_rest_position = points
    fenics_rest_position_rearranged = np.empty(fenics_rest_position.shape)
    fenics_position = np.empty(fenics_rest_position.shape)
    # TODO FEniCS and SOFA do not have the same dofs numbering
    from scipy import spatial

    tree = spatial.KDTree(meshio.read("../meshes/fenics/cube_tetra_1.msh").points)
    for positions, displacement in zip(fenics_rest_position, point_data["f"]):
        fenics_rest_position_rearranged[tree.query(positions)[1]] = positions
        fenics_position[tree.query(positions)[1]] = positions + displacement
    print("initial error:",
          np.mean(fenics_rest_position_rearranged - meshio.read("../meshes/fenics/cube_tetra_1.msh").points))
elif ELEMENT_TYPE == "Hexahedron" and ELEMENT_APPROXIMATION_DEGREE == 1:
    element_sofa = "Hexahedron"
    element_sonics = element_sofa + "_FEniCS"
    element_surface = "Quad"
    mesh = meshio.read("../meshes/fenics/cube_hexa_1.msh")
    indices_sofa = mesh.cells_dict['hexahedron']
    indices_sonics = indices_sofa[:, [0, 1, 3, 2, 4, 5, 7, 6]]
    a = [i.size for i in mesh.cell_data["gmsh:physical"]]
    indices = np.cumsum(a)
    left = mesh.cells_dict["quad"][0:indices[0]]
    right = mesh.cells_dict["quad"][indices[4]:indices[5]]

    with meshio.xdmf.TimeSeriesReader("../meshes/fenics/cube_hexa_1_solution.xdmf") as reader:
        points, cells = reader.read_points_cells()
        t, point_data, cell_data = reader.read_data(0)

    fenics_rest_position = points
    fenics_rest_position_rearranged = np.empty(fenics_rest_position.shape)
    fenics_position = np.empty(fenics_rest_position.shape)
    # TODO FEniCS and SOFA do not have the same dofs numbering
    from scipy import spatial

    tree = spatial.KDTree(mesh.points)
    for positions, displacement in zip(fenics_rest_position, point_data["f"]):
        fenics_rest_position_rearranged[tree.query(positions)[1]] = positions
        fenics_position[tree.query(positions)[1]] = positions + displacement
    print("initial error:", np.mean(fenics_rest_position_rearranged - mesh.points))
elif ELEMENT_TYPE == "Hexahedron" and ELEMENT_APPROXIMATION_DEGREE == 2:
    element_sofa = "Hexahedron20"
    element_sonics = "Hexahedron_FEniCS20"
    element_surface = "Quad8"
    mesh = meshio.read("../meshes/fenics/cube_hexa_2.msh")
    indices_sofa = mesh.cells_dict['hexahedron20']
    indices_sonics = indices_sofa[:, [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]]
    a = [i.size for i in mesh.cell_data["gmsh:physical"]]
    indices = np.cumsum(a)
    left = mesh.cells_dict["quad8"][0:indices[0]]
    right = mesh.cells_dict["quad8"][indices[4]:indices[5]]

    with meshio.xdmf.TimeSeriesReader("../meshes/fenics/cube_hexa_2_solution.xdmf") as reader:
        points, cells = reader.read_points_cells()
        t, point_data, cell_data = reader.read_data(0)

    fenics_rest_position = points
    fenics_rest_position_rearranged = np.empty(fenics_rest_position.shape)
    fenics_position = np.empty(fenics_rest_position.shape)
    # TODO FEniCS and SOFA do not have the same dofs numbering
    from scipy import spatial

    tree = spatial.KDTree(meshio.read("../meshes/fenics/cube_hexa_1.msh").points)
    for positions, displacement in zip(fenics_rest_position, point_data["f"]):
        fenics_rest_position_rearranged[tree.query(positions)[1]] = positions
        fenics_position[tree.query(positions)[1]] = positions + displacement
    print("initial error:",
          np.mean(fenics_rest_position_rearranged - meshio.read("../meshes/fenics/cube_hexa_1.msh").points))


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
                       pluginName="SofaExporter SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")
        root.gravity = [0, 0, 0]

        sonics_node = root.addChild("sonics_node")
        sonics_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-30",
                              absolute_correction_tolerance_threshold="1e-30",
                              absolute_residual_tolerance_threshold="1e-5",
                              relative_residual_tolerance_threshold="1e-10",
                              printLog="1")
        sonics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        self.sonics_mo = sonics_node.addObject('MechanicalObject', name="mo", position=mesh.points.tolist(),
                                               showObject=False, showObjectScale=10, showColor=[0, 1, 0, 1])
        self.sonics_rest_position = np.array(self.sonics_mo.position.value.copy().tolist())

        sonics_node.addObject('CaribouTopology', name='topology', template=element_sonics,
                              indices=indices_sonics.tolist())
        sonics_node.addObject('FixedConstraint', indices=left.tolist())
        sonics_node.addObject("CaribouMass")
        sonics_node.addObject('FEniCS_Material', template=element_sonics, young_modulus="3000",
                              poisson_ratio="0.3", material_name=material, path="/home/..")
        sonics_node.addObject('HyperelasticForcefield_FEniCS', printLog=True)

        fenics_node = root.addChild("fenics_node")
        self.fenics_mo = fenics_node.addObject('MechanicalObject', name="mo", position=fenics_position.tolist(),
                                               showObject=True, showObjectScale=10)

        sofa_node = root.addChild("sofa_node")
        sofa_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-30",
                            absolute_correction_tolerance_threshold="1e-30",
                            absolute_residual_tolerance_threshold="1e-5", relative_residual_tolerance_threshold="1e-10",
                            printLog="1")
        sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        self.sofa_mo = sofa_node.addObject('MechanicalObject', name="mo", position=mesh.points.tolist())
        self.sofa_rest_position = np.array(self.sofa_mo.position.value.copy().tolist())
        sofa_node.addObject('CaribouTopology', name='topology', template=element_sofa,
                            indices=indices_sofa.tolist())
        sofa_node.addObject('FixedConstraint', indices=left.tolist())
        sofa_node.addObject("CaribouMass")
        sofa_node.addObject('CaribouTopology', name='surface', template=element_surface,
                            indices=right.tolist())
        ## For P2 (debug purpose)
        # sofa_node.addObject('CaribouTopology', name='surface', template=element_surface,
        #                     indices=[[4, 20, 97, 21, 101, 102], [29, 4, 97, 31, 102, 103], [20, 5, 98, 22, 104, 105],
        #                              [6, 26, 99, 27, 109, 107],
        #                              [7, 29, 100, 30, 112, 110], [97, 98, 100, 113, 117, 116],
        #                              [97, 20, 98, 101, 105, 113], [98, 23, 99, 106, 108, 114],
        #                              [98, 99, 100, 114, 115, 117], [5, 23, 98, 24, 106, 104], [23, 6, 99, 25, 107, 108],
        #                              [26, 7, 100, 28, 110, 111],
        #                              [99, 26, 100, 109, 111, 115], [29, 97, 100, 103, 116, 112]])
        ## For P1 (debug purpose)
        # sofa_node.addObject('CaribouTopology', name='surface', template=element_surface,
        # indices = [
        #            [4, 40, 12],
        #            [40, 15, 43], [14, 42, 43], [5, 41, 13], [13, 6, 42], [7, 14, 43],
        #            [40, 41, 43], [12, 5, 41], [42, 6, 14], [7, 15, 43], [15, 4, 40],
        #            [41, 42, 43], [40, 12, 41], [42, 41, 13]
        #            ]
        # )

        sofa_node.addObject('TractionForcefield', traction=TRACTION, slope=0, topology='@surface',
                            printLog=True)
        sofa_node.addObject(material + "Material", young_modulus="3000", poisson_ratio="0.3")
        sofa_node.addObject('HyperelasticForcefield', printLog=True)

        return root

    def onSimulationInitDoneEvent(self, event):
        pass

    def onAnimateBeginEvent(self, event):
        pass

    def onAnimateEndEvent(self, event):
        sonics_current_positions = np.array(self.sonics_mo.position.value.copy().tolist())
        sofa_current_positions = np.array(self.sofa_mo.position.value.copy().tolist())

        if ELEMENT_TYPE == "Hexahedron" and ELEMENT_APPROXIMATION_DEGREE == 2:
            sofa_current_positions = sofa_current_positions[np.unique(indices_sofa[:, :8])]
            sonics_current_positions = sonics_current_positions[np.unique(indices_sonics[:, :8])]
        elif ELEMENT_TYPE == "Tetrahedron" and ELEMENT_APPROXIMATION_DEGREE == 2:
            sonics_current_positions = sonics_current_positions[np.unique(indices_sofa[:, :4])]
            sofa_current_positions = sofa_current_positions[np.unique(indices_sonics[:, :4])]

        if self.fenics_mo.position.value.any():
            fenics_current_positions = np.array(self.fenics_mo.position.value.copy().tolist())
            errors = []
            for fenics_current_point, sonics_current_point, sonics_initial_point in zip(fenics_current_positions,
                                                                                        sonics_current_positions,
                                                                                        self.sonics_rest_position):
                if np.linalg.norm(sonics_current_point - sonics_initial_point) != 0:
                    errors.append(np.linalg.norm(fenics_current_point - sonics_current_point) / np.linalg.norm(
                        fenics_current_point - sonics_initial_point))
            mean_error = np.mean(np.array(errors))
            print(f"Relative Mean Error sonics-fenics: {100 * mean_error} %")

        if self.sofa_mo.position.value.any():
            errors = []
            for sofa_current_point, fenics_current_point, sofa_initial_point in zip(sofa_current_positions,
                                                                                    fenics_current_positions,
                                                                                    self.sofa_rest_position):
                if np.linalg.norm(sofa_current_point - sofa_initial_point) != 0:
                    errors.append(np.linalg.norm(sofa_current_point - fenics_current_point) / np.linalg.norm(
                        sofa_current_point - sofa_initial_point))

            mean_error = np.mean(np.array(errors))
            print(f"Relative Mean Error sofa-fenics: {100 * mean_error} %")

            errors = []
            for sofa_current_point, sonics_current_point, sofa_initial_point in zip(sofa_current_positions,
                                                                                    sonics_current_positions,
                                                                                    self.sofa_rest_position):
                if np.linalg.norm(sofa_current_point - sofa_initial_point) != 0:
                    errors.append(np.linalg.norm(sofa_current_point - sonics_current_point) / np.linalg.norm(
                        sofa_current_point - sofa_initial_point))

            print(f"Error sonics-sofa: {sum(errors)}")


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
