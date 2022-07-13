import sys
from pathlib import Path

import Sofa
import SofaCaribou
import numpy as np
from scipy.sparse import csr_matrix, linalg

current_dir = Path(__file__).parent
site_packages_dir = (current_dir / '..' / '..' / 'lib' / 'python3' / 'site-packages').resolve()
sys.path.insert(0, str(site_packages_dir))
print(f'Adding {site_packages_dir} to sys.path')

# TODO: pass the tests for quadratic hexahedron or lower the tolerance
# TODO: better way of doing it ?
TEST_CASES = [
    # ("Kirchhoff_tetra", "SaintVenantKirchhoff", "Tetrahedron"),
    # ("Kirchhoff_tetra10", "SaintVenantKirchhoff", "Tetrahedron10"),
    # ("Kirchhoff_hexa", "SaintVenantKirchhoff", "Hexahedron"),
    # ("Kirchhoff_hexa20", "SaintVenantKirchhoff", "Hexahedron20"),
    # ("NeoHooke_tetra", "NeoHookean", "Tetrahedron"),
    # ("NeoHooke_tetra10", "NeoHookean", "Tetrahedron10"),
    # ("NeoHooke_hexa", "NeoHookean", "Hexahedron"),
    ("NeoHooke_hexa20", "NeoHookean", "Hexahedron20"),
]


def generate_geometry(element):
    hexa_p_0 = [0, 0, 0]
    hexa_p_1 = [1, 0, 0]
    hexa_p_2 = [1, 1, 0]
    hexa_p_3 = [0, 1, 0]
    hexa_p_4 = [0, 0, 1]
    hexa_p_5 = [1, 0, 1]
    hexa_p_6 = [1, 1, 1]
    hexa_p_7 = [0, 1, 1]
    hexa_p_8 = [0.5, 0, 0]
    hexa_p_9 = [1, 0.5, 0]
    hexa_p_10 = [0.5, 1, 0]
    hexa_p_11 = [0, 0.5, 0]
    hexa_p_12 = [0.5, 0, 1]
    hexa_p_13 = [1, 0.5, 1]
    hexa_p_14 = [0.5, 1, 1]
    hexa_p_15 = [0, 0.5, 1]
    hexa_p_16 = [0, 0, 0.5]
    hexa_p_17 = [1, 0, 0.5]
    hexa_p_18 = [1, 1, 0.5]
    hexa_p_19 = [0, 1, 0.5]

    tetra_p_0 = [0, 0, 0]
    tetra_p_1 = [1, 0, 0]
    tetra_p_2 = [0, 1, 0]
    tetra_p_3 = [0, 0, 1]
    tetra_p_4 = [0.5, 0, 0]
    tetra_p_5 = [0.5, 0.5, 0]
    tetra_p_6 = [0, 0.5, 0]
    tetra_p_7 = [0, 0, 0.5]
    tetra_p_8 = [0.5, 0, 0.5]
    tetra_p_9 = [0, 0.5, 0.5]

    if element == "Tetrahedron":
        return element, element, np.array([tetra_p_0, tetra_p_1, tetra_p_2, tetra_p_3]), np.arange(4), np.arange(4)
    elif element == "Tetrahedron10":
        return element, element, np.array(
            [tetra_p_0, tetra_p_1, tetra_p_2, tetra_p_3, tetra_p_4, tetra_p_5, tetra_p_6, tetra_p_7, tetra_p_8,
             tetra_p_9]), np.arange(10), np.array([0, 1, 2, 3, 9, 8, 5, 7, 6, 4])
    elif element == "Hexahedron":
        return element, element + "_FEniCS", np.array(
            [hexa_p_0, hexa_p_1, hexa_p_2, hexa_p_3, hexa_p_4, hexa_p_5, hexa_p_6, hexa_p_7]), np.arange(8), np.array(
            [0, 1, 3, 2, 4, 5, 7, 6])
        # return element, element + "_FEniCS", np.array(
        #     [hexa_p_0, hexa_p_1, hexa_p_2, hexa_p_3, hexa_p_4, hexa_p_5, hexa_p_6, hexa_p_7]), np.arange(8), np.array(
        #     [4, 5, 0, 1, 7, 6, 3, 2])
    elif element == "Hexahedron20":
        return element, "Hexahedron_FEniCS20", np.array(
            [hexa_p_0, hexa_p_1, hexa_p_2, hexa_p_3, hexa_p_4, hexa_p_5, hexa_p_6, hexa_p_7, hexa_p_8, hexa_p_9,
             hexa_p_10, hexa_p_11, hexa_p_12, hexa_p_13, hexa_p_14, hexa_p_15, hexa_p_16, hexa_p_17, hexa_p_18,
             hexa_p_19]), np.arange(20), np.array(
        [4, 5, 0, 1, 7, 6, 3, 2, 12, 16, 15, 17, 13, 8, 11, 9, 14, 19, 18, 10]
        # [0, 1, 3, 2, 4, 5, 7, 6, 8, 11, 16, 9, 17, 10, 19, 18, 12, 15, 13, 14]
        )


class ControlFrame(Sofa.Core.Controller):

    def __init__(self, node):
        Sofa.Core.Controller.__init__(self)
        self.root = self.CreateGraph(node)

    def CreateGraph(self, node):
        node.addObject('DefaultVisualManagerLoop')
        node.addObject('DefaultAnimationLoop')
        node.addObject('VisualStyle', displayFlags="showForceFields showBehaviorModels")
        node.addObject('RequiredPlugin',
                       pluginName="SofaOpenglVisual SofaBaseMechanics SofaBaseTopology SofaSparseSolver SofaImplicitOdeSolver SofaTopologyMapping SofaBoundaryCondition SofaEngine")
        element = "Tetrahedron"
        material = "NeoHookean"
        element_sofa, element_fenics, position, indices_sofa, indices_fenics = generate_geometry(element)
        sofa_node = node.addChild("sofa_node")
        sofa_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-30",
                            absolute_correction_tolerance_threshold="1e-30",
                            absolute_residual_tolerance_threshold="1e-5",
                            relative_residual_tolerance_threshold="1e-10",
                            printLog="1")
        sofa_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        # sofa_node.addObject('StaticODESolver', newton_iterations="20",
        #                     # relative_correction_tolerance_threshold="1e-30",
        #                     #   absolute_correction_tolerance_threshold="1e-30",
        #                     #   absolute_residual_tolerance_threshold="1e-5",
        #                     #   relative_residual_tolerance_threshold="1e-10",
        #                       printLog="1")
        # sofa_node.addObject('LUSolver')
        self.sofa_mo = sofa_node.addObject('MechanicalObject', name="mo", position=position.tolist())
        sofa_node.addObject('CaribouTopology', name='topology', template=element_sofa, indices=indices_sofa.tolist())
        # sofa_node.addObject('CaribouTopology', name='surface', template="Quad8", indices=[4, 5, 6, 7, 12, 13, 14, 15])
        # sofa_node.addObject('CaribouTopology', name='surface', template="Quad", indices=[4, 5, 6, 7])
        sofa_node.addObject('CaribouTopology', name='surface', template="Triangle", indices=[1, 2, 3])

        sofa_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
        # sofa_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
        sofa_node.addObject('FixedConstraint', indices=[0])
        # sofa_node.addObject('ConstantForceField', force=[0, -10, 0], indices=indices_sofa.tolist())
        # sofa_node.addObject("TractionForcefield", traction=[0, 0, 1000], slope=1, topology='@surface', printLog=True)
        # sofa_node.addObject("CaribouMass")
        sofa_node.addObject(material + "Material", young_modulus="3000", poisson_ratio="0.3")
        sofa_node.addObject('HyperelasticForcefield', name="ff", printLog=False)

        self.fenics_node = node.addChild("fenics_node")
        # self.fenics_node.activated=False
        # self.fenics_node.addObject('StaticSolver', newton_iterations="25", relative_correction_tolerance_threshold="1e-30",
        #                       absolute_correction_tolerance_threshold="1e-30",
        #                       absolute_residual_tolerance_threshold="1e-5",
        #                       relative_residual_tolerance_threshold="1e-10",
        #                       printLog="1")
        # self.fenics_node.addObject('SparseLDLSolver', template="CompressedRowSparseMatrixMat3x3d")
        self.fenics_node.addObject('StaticODESolver', newton_iterations="25",
                              # relative_correction_tolerance_threshold="1e-30",
                              #   absolute_correction_tolerance_threshold="1e-30",
                              #   absolute_residual_tolerance_threshold="1e-5",
                              #   relative_residual_tolerance_threshold="1e-10",
                              printLog="1")
        self.fenics_node.addObject('LUSolver')
        self.fenics_mo = self.fenics_node.addObject('MechanicalObject', name="mo", position=position.tolist())
        self.fenics_node.addObject('CaribouTopology', name='topology', template=element_fenics,
                              indices=indices_fenics.tolist())
        self.fenics_node.addObject('BoxROI', name="fixed_roi", box="-7.5 -7.5 -0.9 7.5 7.5 0.1")
        self.fenics_node.addObject('FixedConstraint', indices="@fixed_roi.indices")
        # self.fenics_node.addObject('ConstantForceField', force=[0, -10, 0], indices=indices_fenics.tolist())
        # self.fenics_node.addObject('CaribouTopology', name='surface', template="Quad8",
        #                       indices=[4, 5, 7, 6, 16, 18, 19, 17]
        #                       )
        # self.fenics_node.addObject("TractionForcefield", traction=[0, 0, 1000], slope=1, topology='@surface', printLog=True)
        # self.fenics_node.addObject("CaribouMass")
        self.fenics_node.addObject('FEniCS_Material', template=element_fenics, young_modulus="3000",
                              poisson_ratio="0.3", material_name=material)
        self.fenics_node.addObject('HyperelasticForcefield_FEniCS', name="ff", printLog=False)

        return node

    def onAnimateEndEvent(self, event):
        self.root.gravity = [0, 0, 0]
        fenics_current_positions = np.array(self.fenics_mo.position.value.copy().tolist())
        sofa_current_positions = np.array(self.sofa_mo.position.value.copy().tolist())
        print(fenics_current_positions - sofa_current_positions)


def createScene(node):
    node.addObject(ControlFrame(node))


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
    # K_fenics = csr_matrix(root.fenics_node.ff.K(), copy=True).todense()
    # from scipy import linalg
    # b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000, 0, 0, 1000, 0, 0, 1000, 0, 0, 1000, 0, 0, 0, 0, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 0, 0, 0, 0, 1000, 0, 0, 1000, 0, 0, 1000, 0, 0, 1000, 0, 0, 1000, 0, 0, 1000, 0, 0, 1000])
    #
    # print(b)
    # x = linalg.solve(K_fenics, b.transpose())
    # print(np.dot(K_fenics, np.array(x)))
    # K_sofa = csr_matrix(root.sofa_node.ff.K(), copy=True)
    # print(K_fenics.toarray()[0])
    # print(K_sofa.toarray()[0])
    # print(linalg.norm(K_fenics - K_sofa))

    if not USE_GUI:
        for iteration in range(10):
            Sofa.Simulation.animate(root, root.dt.value)
    else:
        Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
        Sofa.Gui.GUIManager.createGUI(root, __file__)
        Sofa.Gui.GUIManager.SetDimension(1080, 1080)
        Sofa.Gui.GUIManager.MainLoop(root)
        Sofa.Gui.GUIManager.closeGUI()


if __name__ == '__main__':
    main()
